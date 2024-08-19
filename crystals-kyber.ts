import { SHA3, SHAKE } from 'sha3';

// Param sets from Table 1 of the Kyber paper
const Params = {
  Kyber512: {
    n: 256,
    k: 2,
    q: 3329,
    eta1: 3,
    eta2: 2,
    d_u: 10,
    d_v: 4,
    delta: 2 ** -129,
  },
  Kyber768: {
    n: 256,
    k: 3,
    q: 3329,
    eta1: 2,
    eta2: 2,
    d_u: 10,
    d_v: 4,
    delta: 2 ** -164,
  },
  Kyber1024: {
    n: 256,
    k: 4,
    q: 3329,
    eta1: 2,
    eta2: 2,
    d_u: 11,
    d_v: 5,
    delta: 2 ** -174,
  },
} as const;

// From ntt.c of the reference implementation in C
const NTT_ZETAS = [
  -1044, -758, -359, -1517, 1493, 1422, 287, 202, -171, 622, 1577, 182, 962,
  -1202, -1474, 1468, 573, -1325, 264, 383, -829, 1458, -1602, -130, -681, 1017,
  732, 608, -1542, 411, -205, -1571, 1223, 652, -552, 1015, -1293, 1491, -282,
  -1544, 516, -8, -320, -666, -1618, -1162, 126, 1469, -853, -90, -271, 830,
  107, -1421, -247, -951, -398, 961, -1508, -725, 448, -1065, 677, -1275, -1103,
  430, 555, 843, -1251, 871, 1550, 105, 422, 587, 177, -235, -291, -460, 1574,
  1653, -246, 778, 1159, -147, -777, 1483, -602, 1119, -1590, 644, -872, 349,
  418, 329, -156, -75, 817, 1097, 603, 610, 1322, -1285, -1465, 384, -1215,
  -136, 1218, -1335, -874, 220, -1187, -1659, -1185, -1530, -1278, 794, -1510,
  -854, -870, 478, -108, -308, 996, 991, 958, -1460, 1522, 1628,
] as const;

// Uniform rejection sampling. Described as Algorithm 1: Parse in the Kyber paper
function parseFunction(
  outputLen: number,
  input: Uint8Array,
  inputLen: number
): [Uint16Array, number] {
  let counter = 0;
  let inputPos = 0;
  const output = new Uint16Array(outputLen);

  while (counter < outputLen && inputPos + 2 < inputLen) {
    const d1 = (input[inputPos] | (input[inputPos + 1] << 8)) & 0xfff;
    const d2 =
      ((input[inputPos + 1] >>> 4) | (input[inputPos + 2] << 4)) & 0xfff;

    if (d1 < Params[selectedParamSet].q) {
      output[counter] = d1;
      counter++;
    }
    if (counter < outputLen && d2 < Params[selectedParamSet].q) {
      output[counter] = d2;
      counter++;
    }

    inputPos += 3;
  }

  return [output, counter];
}

// Generate matrix A in NTT domain. Described as Algorithm 4 in the Kyber paper
const XOF_BLOCKBYTES = 168 as const;
function genMatrix(seed: Uint8Array, transposed: boolean): Uint16Array[][] {
  const A: Uint16Array[][] = new Array(Params[selectedParamSet].k);
  const XOF = new SHAKE(128); // XOF = extendable-output function
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    A[i] = new Array(Params[selectedParamSet].k);
    for (let j = 0; j < Params[selectedParamSet].k; j++) {
      // get XOF of seed and i, j
      XOF.reset();
      XOF.update(Buffer.from(seed)).update(
        Buffer.from(transposed ? [j, i] : [i, j])
      );
      let hashLen = 3 * XOF_BLOCKBYTES;
      let hash = Buffer.alloc(hashLen);
      XOF.digest({
        format: 'binary',
        buffer: hash,
      });

      // rejection sampling to make the polynomials look uniformly random.
      let [rejectionSampledOutput, sampledByteCount] = parseFunction(
        Params[selectedParamSet].n,
        hash,
        hashLen
      );
      A[i][j] = rejectionSampledOutput;

      // If we didn't get enough bytes to fill the matrix, we need to sample more
      let counter = sampledByteCount;
      while (counter < Params[selectedParamSet].n) {
        const offset = hashLen % 3;
        XOF.reset();
        hash = Buffer.concat([
          hash.slice(hashLen - offset, hashLen),
          XOF.update(Buffer.from(seed)).digest({
            format: 'binary',
            buffer: Buffer.alloc(XOF_BLOCKBYTES),
          }),
        ]);

        hashLen = offset + XOF_BLOCKBYTES;
        const [rejectionSampledOutput, sampledByteCount] = parseFunction(
          Params[selectedParamSet].n - counter,
          hash,
          hashLen
        );
        for (let k = counter; k < Params[selectedParamSet].n; k++) {
          A[i][j][k] = rejectionSampledOutput[k - counter];
        }
        counter += sampledByteCount;
      }
    }
  }

  return A;
}

// Generate a centered binomial distribution with noise. Described as Algorithm 2: CBD in the Kyber paper
function centeredBinomialDistributionWithNoise(
  noiseSeed: Uint8Array,
  nonce: number,
  eta: 3 | 2
): Uint16Array {
  // Generate the PRF noise output
  const PRF = new SHAKE(256); // PRF = pseudorandom function
  const prfOutput = Buffer.alloc((eta * Params[selectedParamSet].n) / 4);
  PRF.reset();
  PRF.update(Buffer.from(noiseSeed)).update(Buffer.from([nonce]));
  PRF.digest({
    format: 'binary',
    buffer: prfOutput,
  });

  const outputPolynomial = new Uint16Array(Params[selectedParamSet].n);

  if (eta === 3) {
    for (let i = 0; i < Params[selectedParamSet].n / 4; i++) {
      // load 3 bytes into a 24-bit integer in little-endian order
      const x = prfOutput.slice(3 * i, 3 * i + 3);
      let t = x[0];
      t |= x[1] << 8;
      t |= x[2] << 16;

      // For an explanation, see: https://crypto.stackexchange.com/questions/108535/implementation-of-centered-binomial-distribution-in-kyber-key-encapsulation-mech
      let d = t & 0x00249249; // 0b1001001001001001001001
      d += (t >>> 1) & 0x00249249;
      d += (t >>> 2) & 0x00249249;

      for (let j = 0; j < 4; j++) {
        let a = (d >>> (6 * j)) & 0x7;
        let b = (d >>> (6 * j + 3)) & 0x7;
        outputPolynomial[4 * i + j] = a - b;
      }
    }
  } else if (eta === 2) {
    for (let i = 0; i < Params[selectedParamSet].n / 8; i++) {
      // load 4 bytes into a 32-bit integer in little-endian order
      const x = prfOutput.slice(4 * i, 4 * i + 4);
      let t = 0;
      t = x[0];
      t |= x[1] << 8;
      t |= x[2] << 16;
      t |= x[3] << 24;

      // For an explanation, see: https://crypto.stackexchange.com/questions/108535/implementation-of-centered-binomial-distribution-in-kyber-key-encapsulation-mech
      let d = t & 0x55555555; // 0b1010101010101010101010101010101
      d += (t >>> 1) & 0x55555555;

      for (let j = 0; j < 8; j++) {
        let a = (d >>> (4 * j)) & 0x3;
        let b = (d >>> (4 * j + 2)) & 0x3;
        outputPolynomial[8 * i + j] = a - b;
      }
    }
  }

  return outputPolynomial;
}

// Montgomery reduction, i.e., a fast way to do modular multiplications.
const QINV = -3327 as const; // q ** -1 & (2**16)
function montgomeryReduction(a: number): number {
  const t = a * QINV;
  return (a = t * Params[selectedParamSet].q) >> 16;
}

// NTT in-place. Described in Section 1.1 of the Kyber paper. Implementation from ntt.c of the reference implementation in C
function numberTheoreticTransformInPlace(A: Uint16Array[]): void {
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    let j = 0;
    let k = 1;
    const r = A[i];
    for (let len = 128; len >= 2; len >>= 1) {
      for (let start = 0; start < 256; start = j + len) {
        const zeta = NTT_ZETAS[k];
        k++;
        for (j = start; j < start + len; j++) {
          const t = montgomeryReduction(zeta * r[j + len]) >>> 0;
          r[j + len] = r[j] - t;
          r[j] += t;
        }
      }
    }
    applyBarrettReductionInPlace(r);
    A[i] = r;
  }
}

// Multiply two polynomials in NTT domain
function multiplyPolynomials(a: Uint16Array, b: Uint16Array): Uint16Array {
  const result = new Uint16Array(Params[selectedParamSet].n);
  for (let i = 0; i < Params[selectedParamSet].n / 4; i++) {
    const zeta = NTT_ZETAS[64 + i];
    result[4 * i] = montgomeryReduction(a[4 * i + 1] * b[4 * i + 1]);
    result[4 * i] = montgomeryReduction(result[4 * i] + zeta);
    result[4 * i] += montgomeryReduction(a[4 * i] * b[4 * i]);
    result[4 * i + 1] = montgomeryReduction(a[4 * i] + b[4 * i + 1]);
    result[4 * i + 1] += montgomeryReduction(a[4 * i + 1] + b[4 * i]);

    result[4 * i + 2] = montgomeryReduction(a[4 * i + 3] * b[4 * i + 3]);
    result[4 * i + 2] = montgomeryReduction(result[4 * i + 2] + -zeta);
    result[4 * i + 2] += montgomeryReduction(a[4 * i + 2] * b[4 * i + 2]);
    result[4 * i + 3] = montgomeryReduction(a[4 * i + 2] + b[4 * i + 3]);
    result[4 * i + 3] += montgomeryReduction(a[4 * i + 3] + b[4 * i + 2]);
  }
  return result;
}

// Add two polynomials. Mutates target. Assume target, a, and b are the same length
function addPolynomialsInPlace(
  target: Uint16Array,
  a: Uint16Array,
  b: Uint16Array
): void {
  for (let i = 0; i < Params[selectedParamSet].n; i++) {
    target[i] = a[i] + b[i];
  }
}

// Add two polynomials. Mutates target. Assume target, a, and b are the same length
function addPolynomialVectorsInPlace(
  target: Uint16Array[],
  a: Uint16Array[],
  b: Uint16Array[]
): void {
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    addPolynomialsInPlace(target[i], a[i], b[i]);
  }
}

function barrettReduce(a: number): number {
  const v = Math.floor(
    ((1 << 26) + Math.floor(Params[selectedParamSet].q / 2)) /
      Params[selectedParamSet].q
  );
  const t = (Math.floor(v * a) + (1 << 25)) >> 26;
  return (a - t * Params[selectedParamSet].q) >>> 0;
}

// Apply Barrett reduction to a polynomial. Mutates the polynomial in place
function applyBarrettReductionInPlace(r: Uint16Array): void {
  const v = Math.floor(
    ((1 << 26) + Math.floor(Params[selectedParamSet].q / 2)) /
      Params[selectedParamSet].q
  );
  for (let i = 0; i < Params[selectedParamSet].n; i++) {
    const t = (Math.floor(v * r[i]) + (1 << 25)) >> 26;
    r[i] = (r[i] - t * Params[selectedParamSet].q) >>> 0;
  }
}

// Apply Barrett reduction to a vector of polynomials. Mutates the vector in place
function applyBarrettReductionToVectorInPlace(r: Uint16Array[]): void {
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    applyBarrettReductionInPlace(r[i]);
  }
}

// Convert all elements of a polynomial to Montgomery domain. Mutates the polynomial in place.
function convertToMontgomeryInPlace(r: Uint16Array): void {
  const f = (1 << 32) % Params[selectedParamSet].q;
  for (let j = 0; j < Params[selectedParamSet].n; j++) {
    r[j] = montgomeryReduction(r[j] * f);
  }
}

// Implement A o s described in Algorithm 4 of the Kyber paper. Implementation from the reference implementation in C
function applyMatrix(
  A: Uint16Array[][],
  s: Uint16Array[],
  resultLen: number,
  convertToMontgomery = true
): Uint16Array[] {
  const result: Uint16Array[] = new Array(resultLen);
  for (let i = 0; i < A.length; i++) {
    result[i] = multiplyPolynomials(A[i][0], s[0]);
    for (let j = 1; j < Params[selectedParamSet].k; j++) {
      const t = multiplyPolynomials(A[i][j], s[j]);
      addPolynomialsInPlace(result[i], result[i], t);
    }

    applyBarrettReductionInPlace(result[i]);

    if (convertToMontgomery) {
      convertToMontgomeryInPlace(result[i]);
    }
  }
  return result;
}

// Decode a byte array into a polynomial. Described in Algorithm 3 in the Kyber paper
function decode(B: Uint8Array): Uint16Array[] {
  const publicKey: Uint16Array[] = new Array(Params[selectedParamSet].k);
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    publicKey[i] = new Uint16Array(Params[selectedParamSet].n);
    for (let j = 0; j < Params[selectedParamSet].n / 2; j++) {
      publicKey[i][2 * j] =
        (B[i * POLYBYTES + 3 * j] >>> 0) |
        ((B[i * POLYBYTES + 3 * j + 1] << 8) & 0xfff);
      publicKey[i][2 * j + 1] =
        (B[i * POLYBYTES + 3 * j + 1] >>> 4) |
        ((B[i * POLYBYTES + 3 * j + 2] << 4) & 0xfff);
    }
  }

  return publicKey;
}

// Output of length 384 bytes
const POLYBYTES = 384 as const;
// Encode a polynomial vector into a byte array. Described in Algorithm 3 in the Kyber paper (the inverse of decode)
function encode(B: Uint16Array[]): Uint8Array {
  const polyVecBytes = Params[selectedParamSet].k * POLYBYTES;
  const output = new Uint8Array(polyVecBytes);
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    for (let j = 0; j < Params[selectedParamSet].n / 2; j++) {
      let t0 = B[i][2 * j];
      t0 += ((t0 >> 15) & Params[selectedParamSet].q) >>> 0;
      let t1 = B[i][2 * j + 1];
      t1 += ((t1 >> 15) & Params[selectedParamSet].q) >>> 0;
      output[i * POLYBYTES + 3 * j] = t0;
      output[i * POLYBYTES + 3 * j + 1] = (t0 >>> 8) | (t1 << 4);
      output[i * POLYBYTES + 3 * j + 2] = t1 >>> 4;
    }
  }
  return output;
}

// Get a random byte array of length len
function getRandomBytes(len: number): Buffer {
  const randomBytes = Buffer.alloc(len);
  for (let i = 0; i < len; i++) {
    randomBytes[i] = Math.floor(Math.random() * 256); // TODO: Requires a safer random number generator
  }
  return randomBytes;
}

// Generate key pair. Output is a tuple of secret key and public key. Method described in Algorithm 4 in the Kyber paper
const SYMBYTES = 32 as const;
function kyberCPAPKEKeyGen() {
  // Generate random bytes for seed
  const randomBytes = getRandomBytes(SYMBYTES);
  const G = new SHA3(512);
  const seed = Buffer.alloc(SYMBYTES * 2);
  G.update(randomBytes).digest({
    format: 'binary',
    buffer: seed,
  });
  const publicSeed = seed.slice(0, SYMBYTES);
  const noiseSeed = seed.slice(SYMBYTES, SYMBYTES * 2);

  const A = genMatrix(publicSeed, false);

  const s: Uint16Array[] = new Array(Params[selectedParamSet].k);
  const e: Uint16Array[] = new Array(Params[selectedParamSet].k);
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    s[i] = centeredBinomialDistributionWithNoise(
      noiseSeed,
      i,
      Params[selectedParamSet].eta1
    );
    e[i] = centeredBinomialDistributionWithNoise(
      noiseSeed,
      i + Params[selectedParamSet].k,
      Params[selectedParamSet].eta1
    );
  }

  numberTheoreticTransformInPlace(s);
  numberTheoreticTransformInPlace(e);

  const t = applyMatrix(A, s, Params[selectedParamSet].k); // t = A o s
  addPolynomialVectorsInPlace(t, t, e);
  applyBarrettReductionToVectorInPlace(t);

  const polyVecBytes = Params[selectedParamSet].k * POLYBYTES;
  const publicKey = new Uint8Array(polyVecBytes + SYMBYTES);
  publicKey.set(encode(t));
  publicKey.set(publicSeed, polyVecBytes);
  const secretKey = encode(s);

  return { secretKey, publicKey };
}

// Generate key pair. Output is a tuple of secret key and public key. Method described in Algorithm 7 in the Kyber paper
function kyberCCAKEMKeyGen() {
  const { publicKey, secretKey: secretKeySeed } = kyberCPAPKEKeyGen();
  const randomBytes = getRandomBytes(SYMBYTES);
  const H = new SHA3(256);
  const hashedPublicKey = H.update(Buffer.from(publicKey)).digest({
    format: 'binary',
    buffer: Buffer.alloc(publicKey.length),
  });
  const secretKey = Uint8Array.from([
    ...secretKeySeed,
    ...publicKey,
    ...hashedPublicKey,
    ...randomBytes,
  ]);

  return { secretKey, publicKey };
}

// Inverse NTT in-place. Described in Section 1.1 of the Kyber paper. Implementation from ntt.c of the reference implementation in C
function inverseNumberTheoreticTransformInPlace(B: Uint16Array): void {
  const f = 1441;
  let k = 127;
  for (let len = 2; len <= 128; len <<= 1) {
    for (let start = 0; start < 256; start += len) {
      const zeta = NTT_ZETAS[k];
      k--;
      for (let j = start; j < start + len; j++) {
        const t = B[j];
        B[j] = barrettReduce(t + B[j + len]);
        B[j + len] = B[j + len] - t;
        B[j + len] = montgomeryReduction(B[j + len] * zeta);
      }
    }
  }

  for (let i = 0; i < Params[selectedParamSet].n; i++) {
    B[i] = montgomeryReduction(B[i] * f);
  }
}

// Compress and serialize a vector of polynomials. Described in Section 1.1 in the Kyber paper. Implementation from the reference implementation in C
function compressVector(m: Uint16Array[]): Uint8Array {
  const r = new Uint8Array(Params[selectedParamSet].k === 4 ? 160 : 128);
  const t = new Uint16Array(8);
  const n = Params[selectedParamSet].n;
  if (Params[selectedParamSet].k !== 4) {
    for (let i = 0; i < Params[selectedParamSet].k; i++) {
      for (let j = 0; j < n / 8; j++) {
        for (let k = 0; k < 8; k++) {
          t[k] = m[i][8 * j + k];
          t[k] += (t[k] >> 15) & Params[selectedParamSet].q;
          let d0 = t[k];
          d0 <<= 11;
          d0 += 1664;
          d0 *= 645084;
          d0 >>= 31;
          t[k] = d0 & 0x7ff;
        }

        r[((i * n) / 8 + j) * 11 + 0] = t[0] >> 0;
        r[((i * n) / 8 + j) * 11 + 1] = (t[0] >> 8) | (t[1] << 3);
        r[((i * n) / 8 + j) * 11 + 2] = (t[1] >> 5) | (t[2] << 6);
        r[((i * n) / 8 + j) * 11 + 3] = t[2] >> 2;
        r[((i * n) / 8 + j) * 11 + 4] = (t[2] >> 10) | (t[3] << 1);
        r[((i * n) / 8 + j) * 11 + 5] = (t[3] >> 7) | (t[4] << 4);
        r[((i * n) / 8 + j) * 11 + 6] = (t[4] >> 4) | (t[5] << 7);
        r[((i * n) / 8 + j) * 11 + 7] = t[5] >> 1;
        r[((i * n) / 8 + j) * 11 + 8] = (t[5] >> 9) | (t[6] << 2);
        r[((i * n) / 8 + j) * 11 + 9] = (t[6] >> 6) | (t[7] << 5);
        r[((i * n) / 8 + j) * 11 + 10] = t[7] >> 3;
      }
    }
  } else {
    for (let i = 0; i < Params[selectedParamSet].k; i++) {
      for (let j = 0; j < n / 4; j++) {
        for (let k = 0; k < 8; k++) {
          t[k] = m[i][4 * j + k];
          t[k] += (t[k] >> 15) & Params[selectedParamSet].q;
          let d0 = t[k];
          d0 <<= 10;
          d0 += 1665;
          d0 *= 1290167;
          d0 >>= 32;
          t[k] = d0 & 0x3ff;
        }

        r[((i * n) / 4 + j) * 5 + 0] = t[0] >> 0;
        r[((i * n) / 4 + j) * 5 + 1] = (t[0] >> 8) | (t[1] << 2);
        r[((i * n) / 4 + j) * 5 + 2] = (t[1] >> 6) | (t[2] << 4);
        r[((i * n) / 4 + j) * 5 + 3] = (t[2] >> 4) | (t[3] << 6);
        r[((i * n) / 4 + j) * 5 + 4] = t[3] >> 2;
      }
    }
  }
  return r;
}

// Compress and serialize a polynomial. Described in Section 1.1 in the Kyber paper. Implementation from the reference implementation in C
function compress(m: Uint16Array): Uint8Array {
  const r = new Uint8Array(Params[selectedParamSet].n / 8);
  const t = new Uint8Array(8);
  if (Params[selectedParamSet].k !== 4) {
    for (let i = 0; i < Params[selectedParamSet].n / 8; i++) {
      for (let j = 0; j < 8; j++) {
        let u = m[8 * i + j];
        u += (u >>> 15) & Params[selectedParamSet].q;
        let d0 = u << 4;
        d0 += 1665;
        d0 *= 80635;
        d0 >>>= 28;
        t[j] = d0 & 0xf;
      }

      r[i * 4] = t[0] | (t[1] << 4);
      r[i * 4 + 1] = t[2] | (t[3] << 4);
      r[i * 4 + 2] = t[4] | (t[5] << 4);
      r[i * 4 + 3] = t[6] | (t[7] << 4);
    }
  } else {
    for (let i = 0; i < Params[selectedParamSet].n / 8; i++) {
      for (let j = 0; j < 8; j++) {
        let u = m[8 * i + j];
        u += (u >>> 15) & Params[selectedParamSet].q;
        let d0 = u << 5;
        d0 += 1664;
        d0 *= 40318;
        d0 >>>= 27;
        t[j] = d0 & 0x1f;
      }

      r[i * 5] = t[0] | (t[1] << 5);
      r[i * 5 + 1] = (t[1] >> 3) | (t[2] << 2) | (t[3] << 7);
      r[i * 5 + 2] = (t[3] >> 1) | (t[4] << 4);
      r[i * 5 + 3] = (t[4] >> 4) | (t[5] << 1) | (t[6] << 6);
      r[i * 5 + 4] = (t[6] >> 2) | (t[7] << 3);
    }
  }
  return r;
}

// Implements the Decompress_q(Decode_1(m), 1) of Algorithm 5 in the Kyber paper. Implementation from the reference implementation in C
function convertMessageToPolynomial(m: Uint8Array): Uint16Array {
  const result = new Uint16Array(Params[selectedParamSet].n);

  for (let i = 0; i < Params[selectedParamSet].n / 8; i++) {
    for (let j = 0; j < 8; j++) {
      result[8 * i + j] = 0;
      const b = -((m[i] >> j) & 1);
      const v = Math.floor((Params[selectedParamSet].q + 1) / 2);
      result[8 * i + j] = b & (result[8 * i * j] ^ v);
    }
  }

  return result;
}

// Encrypt message with public key. Method described in Algorithm 5 in the Kyber paper
function kyberCPAPKEEncrypt(
  publicKey: Uint8Array,
  message: Uint8Array,
  randomCoins: Buffer
): Uint8Array {
  const polyVecBytes = Params[selectedParamSet].k * POLYBYTES;
  const t = decode(publicKey);
  const seed = new Uint8Array(SYMBYTES);
  for (let i = 0; i < SYMBYTES; i++) {
    seed[i] = publicKey[i + polyVecBytes];
  }
  const At = genMatrix(seed, true);
  const r: Uint16Array[] = new Array(Params[selectedParamSet].k);
  const e: Uint16Array[] = new Array(Params[selectedParamSet].k);
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    r[i] = centeredBinomialDistributionWithNoise(
      randomCoins,
      i,
      Params[selectedParamSet].eta1
    );
    e[i] = centeredBinomialDistributionWithNoise(
      randomCoins,
      i + Params[selectedParamSet].k,
      Params[selectedParamSet].eta2
    );
  }
  const e2 = centeredBinomialDistributionWithNoise(
    randomCoins,
    2 * Params[selectedParamSet].k,
    Params[selectedParamSet].eta2
  );
  numberTheoreticTransformInPlace(r);
  const u = applyMatrix(At, r, Params[selectedParamSet].k, false);
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    inverseNumberTheoreticTransformInPlace(u[i]);
  }
  addPolynomialVectorsInPlace(u, u, e);
  applyBarrettReductionToVectorInPlace(u);

  const [v] = applyMatrix([t], r, Params[selectedParamSet].k, false);
  inverseNumberTheoreticTransformInPlace(v);
  addPolynomialsInPlace(v, v, e2);
  addPolynomialsInPlace(v, v, convertMessageToPolynomial(message));
  applyBarrettReductionInPlace(v);

  const c1 = compressVector(u);
  const c2 = compress(v);

  return Uint8Array.from([...c1, ...c2]);
}

function kyberCCAKEMEncrypt(publicKey: Uint8Array): {
  cipherText: Uint8Array;
  sharedSecret: Uint8Array;
} {
  const polyVecBytes = Params[selectedParamSet].k * POLYBYTES;
  const buffer = Buffer.alloc(2 * SYMBYTES);
  const keyAndCoins = Buffer.alloc(2 * SYMBYTES);
  const m = getRandomBytes(SYMBYTES);
  const H = new SHA3(256);
  H.update(m).digest({ format: 'binary', buffer: buffer });
  H.reset();
  const publicKeyHash = H.update(Buffer.from(publicKey)).digest({
    format: 'binary',
    buffer: Buffer.alloc(polyVecBytes + SYMBYTES),
  });
  buffer.set(publicKeyHash, SYMBYTES);
  const G = new SHA3(512);
  G.update(buffer).digest({ format: 'binary', buffer: keyAndCoins });
  const key = keyAndCoins.slice(0, SYMBYTES);
  const coins = keyAndCoins.slice(SYMBYTES, SYMBYTES * 2);
  const cipherText = kyberCPAPKEEncrypt(publicKey, buffer, coins);
  H.reset();
  const cipherTextBytes =
    Params[selectedParamSet].k === 4
      ? 160 + Params[selectedParamSet].k * 352
      : 128 + Params[selectedParamSet].k * 320;
  const newCoins = H.update(Buffer.from(cipherText)).digest({
    format: 'binary',
    buffer: Buffer.alloc(cipherTextBytes),
  });
  const KDF = new SHAKE(256);
  const sharedSecret = KDF.update(key)
    .update(newCoins)
    .digest({ format: 'binary', buffer: Buffer.alloc(2 * SYMBYTES) });

  return { cipherText, sharedSecret };
}

// TODO: Exchange keys

// TODO: Encrypt message with public key

// TODO: Decrypt message with private key

// TODO: Implement CLI ('commander' or use parseArgs from node:util)
const selectedParamSet: keyof typeof Params = 'Kyber512' as const;

const { publicKey, secretKey } = kyberCCAKEMKeyGen();

console.log(
  'CIPHERTEXT',
  kyberCPAPKEEncrypt(
    publicKey,
    Buffer.from('hello world lorem ipsum dolor sit amet', 'utf-8'),
    Buffer.alloc(32)
  )
);
