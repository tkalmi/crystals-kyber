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
  },
  Kyber768: {
    n: 256,
    k: 3,
    q: 3329,
    eta1: 2,
    eta2: 2,
    d_u: 10,
    d_v: 4,
  },
  Kyber1024: {
    n: 256,
    k: 4,
    q: 3329,
    eta1: 2,
    eta2: 2,
    d_u: 11,
    d_v: 5,
  },
} as const;

// From the original TS implementation
const NTT_ZETAS = [
  2285, 2571, 2970, 1812, 1493, 1422, 287, 202, 3158, 622, 1577, 182, 962, 2127,
  1855, 1468, 573, 2004, 264, 383, 2500, 1458, 1727, 3199, 2648, 1017, 732, 608,
  1787, 411, 3124, 1758, 1223, 652, 2777, 1015, 2036, 1491, 3047, 1785, 516,
  3321, 3009, 2663, 1711, 2167, 126, 1469, 2476, 3239, 3058, 830, 107, 1908,
  3082, 2378, 2931, 961, 1821, 2604, 448, 2264, 677, 2054, 2226, 430, 555, 843,
  2078, 871, 1550, 105, 422, 587, 177, 3094, 3038, 2869, 1574, 1653, 3083, 778,
  1159, 3182, 2552, 1483, 2727, 1119, 1739, 644, 2457, 349, 418, 329, 3173,
  3254, 817, 1097, 603, 610, 1322, 2044, 1864, 384, 2114, 3193, 1218, 1994,
  2455, 220, 2142, 1670, 2144, 1799, 2051, 794, 1819, 2475, 2459, 478, 3221,
  3021, 996, 991, 958, 1869, 1522, 1628,
] as const;

// From the original TS implementation
const NTT_ZETAS_INV = [
  1701, 1807, 1460, 2371, 2338, 2333, 308, 108, 2851, 870, 854, 1510, 2535,
  1278, 1530, 1185, 1659, 1187, 3109, 874, 1335, 2111, 136, 1215, 2945, 1465,
  1285, 2007, 2719, 2726, 2232, 2512, 75, 156, 3000, 2911, 2980, 872, 2685,
  1590, 2210, 602, 1846, 777, 147, 2170, 2551, 246, 1676, 1755, 460, 291, 235,
  3152, 2742, 2907, 3224, 1779, 2458, 1251, 2486, 2774, 2899, 1103, 1275, 2652,
  1065, 2881, 725, 1508, 2368, 398, 951, 247, 1421, 3222, 2499, 271, 90, 853,
  1860, 3203, 1162, 1618, 666, 320, 8, 2813, 1544, 282, 1838, 1293, 2314, 552,
  2677, 2106, 1571, 205, 2918, 1542, 2721, 2597, 2312, 681, 130, 1602, 1871,
  829, 2946, 3065, 1325, 2756, 1861, 1474, 1202, 2367, 3147, 1752, 2707, 171,
  3127, 3042, 1907, 1836, 1517, 359, 758, 1441,
] as const;

// Converts a number to a 16-bit integer
const int16 = new Int16Array(1);
function toInt16(n: number): number {
  int16[0] = n;
  return int16[0];
}

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
        Buffer.from(transposed ? [i, j] : [j, i])
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
const QINV = 62209 as const; // ((q ** -1 & (2**16)) >>> 0) & 0xffff
function montgomeryReduction(a: number): number {
  const t = toInt16(((a >>> 0) & 0xffff) * QINV);
  return toInt16((a - t * Params[selectedParamSet].q) >> 16);
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
          const rjlenInt = toInt16(r[j + len]);
          const t = montgomeryReduction(zeta * rjlenInt);
          r[j + len] = r[j] - t;
          r[j] += t;
        }
      }
    }
    applyBarrettReductionToPolynomialInPlace(r);
    A[i] = r;
  }
}

// Multiply two polynomials in NTT domain
function multiplyPolynomials(a: Uint16Array, b: Uint16Array): Uint16Array {
  const result = new Uint16Array(Params[selectedParamSet].n);
  for (let i = 0; i < Params[selectedParamSet].n / 4; i++) {
    const zeta = NTT_ZETAS[64 + i];
    const a0 = a[4 * i];
    const b0 = b[4 * i];
    const a1 = a[4 * i + 1];
    const b1 = b[4 * i + 1];
    let tmp = montgomeryReduction(a1 * b1);
    result[4 * i] = montgomeryReduction(tmp * zeta);
    result[4 * i] += montgomeryReduction(a0 * b0);
    result[4 * i + 1] = montgomeryReduction(a0 * b1);
    result[4 * i + 1] += montgomeryReduction(a1 * b0);

    const a2 = a[4 * i + 2];
    const b2 = b[4 * i + 2];
    const a3 = a[4 * i + 3];
    const b3 = b[4 * i + 3];
    tmp = montgomeryReduction(a3 * b3);
    result[4 * i + 2] = montgomeryReduction(tmp * -zeta);
    result[4 * i + 2] += montgomeryReduction(a2 * b2);
    result[4 * i + 3] = montgomeryReduction(a2 * b3);
    result[4 * i + 3] += montgomeryReduction(a3 * b2);
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
  const aInt16 = toInt16(a);
  const v = 20159; // int16(((1<<26) + 3329/2)/3329)
  let t = toInt16((v * aInt16) >> 26);
  t = toInt16(t * Params[selectedParamSet].q);
  return aInt16 - t;
}

// Apply Barrett reduction to a polynomial. Mutates the polynomial in place
function applyBarrettReductionToPolynomialInPlace(r: Uint16Array): void {
  for (let i = 0; i < Params[selectedParamSet].n; i++) {
    r[i] = barrettReduce(r[i]);
  }
}

// Apply Barrett reduction to a vector of polynomials. Mutates the vector in place
function applyBarrettReductionToPolynomialVectorInPlace(
  r: Uint16Array[]
): void {
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    applyBarrettReductionToPolynomialInPlace(r[i]);
  }
}

// Convert all elements of a polynomial to Montgomery domain. Mutates the polynomial in place.
function convertPolynomialToMontgomeryInPlace(r: Uint16Array): void {
  const f = 1353; // (1ULL << 32) % Params[selectedParamSet].q
  for (let j = 0; j < Params[selectedParamSet].n; j++) {
    r[j] = montgomeryReduction(r[j] * f);
  }
}

// Implement A o s described in Algorithm 4 of the Kyber paper. Multiplies elements of A and s in NTT domain, accumulate into R and multiply by 2^-16. Implementation from the reference implementation in C.
function multiplyPolynomialMatrixAndVector(
  A: Uint16Array[][],
  s: Uint16Array[],
  resultLen: number,
  convertToMontgomery: boolean
): Uint16Array[] {
  const result: Uint16Array[] = new Array(resultLen);
  for (let i = 0; i < A.length; i++) {
    result[i] = multiplyPolynomials(A[i][0], s[0]);
    for (let j = 1; j < Params[selectedParamSet].k; j++) {
      const t = multiplyPolynomials(A[i][j], s[j]);
      addPolynomialsInPlace(result[i], result[i], t);
    }

    applyBarrettReductionToPolynomialInPlace(result[i]);

    if (convertToMontgomery) {
      convertPolynomialToMontgomeryInPlace(result[i]);
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
export function kyberCPAPKEKeyGen() {
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

  const t = multiplyPolynomialMatrixAndVector(
    A,
    s,
    Params[selectedParamSet].k,
    true
  ); // t = A o s
  addPolynomialVectorsInPlace(t, t, e);
  applyBarrettReductionToPolynomialVectorInPlace(t);

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

// Inverse NTT in-place. Described in Section 1.1 of the Kyber paper. Implementation from the original TS implementation
function inverseNumberTheoreticTransformInPlace(B: Uint16Array): void {
  const f = 1441;
  let k = 0;
  let j = 0;
  for (let len = 2; len <= 128; len <<= 1) {
    for (let start = 0; start < 256; start = j + len) {
      const zeta = NTT_ZETAS_INV[k];
      k++;
      for (j = start; j < start + len; j++) {
        const t = toInt16(B[j]);
        B[j] = barrettReduce(t + toInt16(B[j + len]));
        const tmp = t - toInt16(B[j + len]);
        B[j + len] = montgomeryReduction(tmp * zeta);
      }
    }
  }

  for (let i = 0; i < 256; i++) {
    B[i] = montgomeryReduction(toInt16(B[i]) * f);
  }
}

// Compress and serialize a vector of polynomials. Described in Section 1.1 in the Kyber paper.
function compressPolynomialVector(m: Uint16Array[]): Uint8Array {
  const r = new Uint8Array(
    Params[selectedParamSet].k === 4
      ? Params[selectedParamSet].k * 352
      : Params[selectedParamSet].k * 320
  );
  const n = Params[selectedParamSet].n;
  let rPos = 0;
  if (Params[selectedParamSet].k === 4) {
    const t = new Uint16Array(8);
    for (let i = 0; i < Params[selectedParamSet].k; i++) {
      for (let j = 0; j < n / 8; j++) {
        for (let k = 0; k < 8; k++) {
          const x = m[i][8 * j + k];
          t[k] =
            (Math.round((2048 / Params[selectedParamSet].q) * x) >>> 0) % 2048;
        }

        r[rPos + 0] = t[0] >> 0;
        r[rPos + 1] = (t[0] >> 8) | (t[1] << 3);
        r[rPos + 2] = (t[1] >> 5) | (t[2] << 6);
        r[rPos + 3] = t[2] >> 2;
        r[rPos + 4] = (t[2] >> 10) | (t[3] << 1);
        r[rPos + 5] = (t[3] >> 7) | (t[4] << 4);
        r[rPos + 6] = (t[4] >> 4) | (t[5] << 7);
        r[rPos + 7] = t[5] >> 1;
        r[rPos + 8] = (t[5] >> 9) | (t[6] << 2);
        r[rPos + 9] = (t[6] >> 6) | (t[7] << 5);
        r[rPos + 10] = t[7] >> 3;

        rPos += 11;
      }
    }
  } else {
    const t = new Uint16Array(4);
    for (let i = 0; i < Params[selectedParamSet].k; i++) {
      for (let j = 0; j < n / 4; j++) {
        for (let k = 0; k < 4; k++) {
          const x = m[i][4 * j + k];
          t[k] =
            (Math.round((1024 / Params[selectedParamSet].q) * x) >>> 0) % 1024;
        }

        r[rPos + 0] = t[0] >> 0;
        r[rPos + 1] = (t[0] >> 8) | (t[1] << 2);
        r[rPos + 2] = (t[1] >> 6) | (t[2] << 4);
        r[rPos + 3] = (t[2] >> 4) | (t[3] << 6);
        r[rPos + 4] = t[3] >> 2;

        rPos += 5;
      }
    }
  }
  return r;
}

// Compress and serialize a polynomial. Described in Section 1.1 in the Kyber paper.
function compressPolynomial(m: Uint16Array): Uint8Array {
  const r = new Uint8Array(Params[selectedParamSet].k === 4 ? 160 : 128);
  const t = new Uint8Array(8);
  if (Params[selectedParamSet].k !== 4) {
    for (let i = 0; i < Params[selectedParamSet].n / 8; i++) {
      for (let j = 0; j < 8; j++) {
        const x = m[8 * i + j];
        t[j] = (Math.round((16 / Params[selectedParamSet].q) * x) >>> 0) % 16;
      }

      r[i * 4] = t[0] | (t[1] << 4);
      r[i * 4 + 1] = t[2] | (t[3] << 4);
      r[i * 4 + 2] = t[4] | (t[5] << 4);
      r[i * 4 + 3] = t[6] | (t[7] << 4);
    }
  } else {
    for (let i = 0; i < Params[selectedParamSet].n / 8; i++) {
      for (let j = 0; j < 8; j++) {
        const x = m[8 * i + j];
        t[j] = (Math.round((32 / Params[selectedParamSet].q) * x) >>> 0) % 32;
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
      const b = -((((m[i] >> j) & 1) >>> 0) & 0xffff);
      const v = Math.floor((Params[selectedParamSet].q + 1) / 2);
      result[8 * i + j] = b & (0 ^ v);
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
  const u = multiplyPolynomialMatrixAndVector(
    At,
    r,
    Params[selectedParamSet].k,
    false
  );
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    inverseNumberTheoreticTransformInPlace(u[i]);
  }
  addPolynomialVectorsInPlace(u, u, e);
  applyBarrettReductionToPolynomialVectorInPlace(u);

  const [v] = multiplyPolynomialMatrixAndVector(
    [t],
    r,
    Params[selectedParamSet].k,
    false
  );
  inverseNumberTheoreticTransformInPlace(v);
  addPolynomialsInPlace(v, v, e2);
  addPolynomialsInPlace(v, v, convertMessageToPolynomial(message));
  applyBarrettReductionToPolynomialInPlace(v);

  const c1 = compressPolynomialVector(u);
  const c2 = compressPolynomial(v);

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
  buffer.set(publicKeyHash.slice(0, SYMBYTES), SYMBYTES);
  const G = new SHA3(512);
  G.update(buffer).digest({ format: 'binary', buffer: keyAndCoins });
  const key = keyAndCoins.slice(0, SYMBYTES);
  const coins = keyAndCoins.slice(SYMBYTES, SYMBYTES * 2);
  const cipherText = kyberCPAPKEEncrypt(
    publicKey,
    buffer.slice(0, SYMBYTES),
    coins
  );
  H.reset();
  const newCoins = H.update(Buffer.from(cipherText)).digest();
  const KDF = new SHAKE(256);
  const sharedSecret = KDF.update(key)
    .update(newCoins)
    .digest({ format: 'binary', buffer: Buffer.alloc(SYMBYTES) });

  return { cipherText, sharedSecret };
}

function decompressPolynomialVector(a: Uint8Array): Uint16Array[] {
  const result = new Array(Params[selectedParamSet].k);
  let aPos = 0;
  if (Params[selectedParamSet].k === 4) {
    const t = new Uint16Array(8);
    for (let i = 0; i < Params[selectedParamSet].k; i++) {
      result[i] = new Uint16Array(Params[selectedParamSet].n);
      for (let j = 0; j < Params[selectedParamSet].n / 8; j++) {
        t[0] = a[aPos + 0] | (a[aPos + 1] << 8);
        t[1] = (a[aPos + 1] >> 3) | (a[aPos + 2] << 5);
        t[2] = (a[aPos + 2] >> 6) | (a[aPos + 3] << 2) | (a[aPos + 4] << 10);
        t[3] = (a[aPos + 4] >> 1) | (a[aPos + 5] << 7);
        t[4] = (a[aPos + 5] >> 4) | (a[aPos + 6] << 4);
        t[5] = (a[aPos + 6] >> 7) | (a[aPos + 7] << 1) | (a[aPos + 8] << 9);
        t[6] = (a[aPos + 8] >> 2) | (a[aPos + 9] << 6);
        t[7] = (a[aPos + 9] >> 5) | (a[aPos + 10] << 3);
        aPos += 11;

        for (let k = 0; k < 8; k++) {
          result[i][8 * j + k] =
            ((t[k] & 0x7ff) * Params[selectedParamSet].q + 1024) >> 11;
        }
      }
    }
  } else {
    const t = new Uint16Array(4);
    for (let i = 0; i < Params[selectedParamSet].k; i++) {
      result[i] = new Uint16Array(Params[selectedParamSet].n);
      for (let j = 0; j < Params[selectedParamSet].n / 4; j++) {
        t[0] = a[aPos + 0] | (a[aPos + 1] << 8);
        t[1] = (a[aPos + 1] >> 2) | (a[aPos + 2] << 6);
        t[2] = (a[aPos + 2] >> 4) | (a[aPos + 3] << 4);
        t[3] = (a[aPos + 3] >> 6) | (a[aPos + 4] << 2);
        aPos += 5;

        for (let k = 0; k < 4; k++) {
          result[i][4 * j + k] =
            ((t[k] & 0x3ff) * Params[selectedParamSet].q + 512) >> 10;
        }
      }
    }
  }

  return result;
}

function decompressPolynomial(a: Uint8Array): Uint16Array {
  const result = new Uint16Array(Params[selectedParamSet].n);

  let aPos = 0;
  if (Params[selectedParamSet].k === 4) {
    const t = new Uint8Array(8);
    for (let i = 0; i < Params[selectedParamSet].n / 8; i++) {
      t[0] = a[aPos];
      t[1] = (a[aPos] >> 5) | (a[aPos + 1] << 3);
      t[2] = a[aPos + 1] >> 2;
      t[3] = (a[aPos + 1] >> 7) | (a[aPos + 2] << 1);
      t[4] = (a[aPos + 2] >> 4) | (a[aPos + 3] << 4);
      t[5] = a[aPos + 3] >> 1;
      t[6] = (a[aPos + 3] >> 6) | (a[aPos + 4] << 2);
      t[7] = a[aPos + 4] >> 3;
      aPos += 5;

      for (let j = 0; j < 8; j++) {
        result[8 * i + j] =
          ((t[j] & 31) * Params[selectedParamSet].q + 16) >> 5;
      }
    }
  } else {
    for (let i = 0; i < Params[selectedParamSet].n / 2; i++) {
      result[2 * i] = ((a[aPos] & 15) * Params[selectedParamSet].q + 8) >> 4;
      result[2 * i + 1] =
        ((a[aPos] >> 4) * Params[selectedParamSet].q + 8) >> 4;
      aPos++;
    }
  }

  return result;
}

function subtractPolynomialsInPlace(
  target: Uint16Array,
  a: Uint16Array,
  b: Uint16Array
): void {
  for (let i = 0; i < Params[selectedParamSet].n; i++) {
    target[i] = a[i] - b[i];
  }
}

function convertPolynomialToMessage(a: Uint16Array): Uint8Array {
  const result = new Uint8Array(SYMBYTES);
  for (let i = 0; i < Params[selectedParamSet].n / 8; i++) {
    result[i] = 0;
    for (let j = 0; j < 8; j++) {
      let t = a[8 * i + j];
      t <<= 1;
      t += 1665;
      t *= 80635;
      t >>= 28;
      t &= 1;
      result[i] |= t << j;
    }
  }

  return result;
}

// Decrypt message with private key
function kyberCPAPKEDecrypt(
  secretKey: Uint8Array,
  cipherText: Uint8Array
): Uint8Array {
  const u = decompressPolynomialVector(cipherText);
  const POLYVECCOMPRESSEDBYTES = Params[selectedParamSet].k === 4 ? 352 : 320;
  const v = decompressPolynomial(
    cipherText.slice(Params[selectedParamSet].k * POLYVECCOMPRESSEDBYTES)
  );
  const s = decode(secretKey);

  numberTheoreticTransformInPlace(u);
  const [mp] = multiplyPolynomialMatrixAndVector(
    [s],
    u,
    Params[selectedParamSet].k,
    false
  );
  inverseNumberTheoreticTransformInPlace(mp);

  subtractPolynomialsInPlace(mp, v, mp);
  applyBarrettReductionToPolynomialInPlace(mp);

  return convertPolynomialToMessage(mp);
}

function kyberCCAKEMDecrypt(
  secretKey: Uint8Array,
  cipherText: Uint8Array
): Uint8Array {
  const privateKeyLength =
    (12 * Params[selectedParamSet].k * Params[selectedParamSet].n) / 8;
  const privateKey = secretKey.slice(0, privateKeyLength);
  const publicKey = secretKey.slice(
    privateKeyLength,
    privateKeyLength * 2 + SYMBYTES
  );
  const h = secretKey.slice(
    2 * privateKeyLength + SYMBYTES,
    2 * privateKeyLength + SYMBYTES * 2
  );
  const z = secretKey.slice(
    2 * privateKeyLength + SYMBYTES * 2,
    2 * privateKeyLength + SYMBYTES * 3
  );

  const m = kyberCPAPKEDecrypt(privateKey, cipherText);
  const G = new SHA3(512);
  const keyAndCoins = G.update(Buffer.from(m)).update(Buffer.from(h)).digest();
  const key = keyAndCoins.slice(0, SYMBYTES);
  const coins = keyAndCoins.slice(SYMBYTES);
  const comparisonCipher = kyberCPAPKEEncrypt(publicKey, Buffer.from(m), coins);
  const KDF = new SHAKE(256);
  const H = new SHA3(256);
  const cipherHash = H.update(Buffer.from(comparisonCipher)).digest();
  if (comparisonCipher.toString() === cipherText.toString()) {
    return KDF.update(key).update(cipherHash).digest();
  } else {
    return KDF.update(Buffer.from(z)).update(cipherHash).digest();
  }
}

// TODO: Implement CLI ('commander' or use parseArgs from node:util)
// This is a KEM algorithm. So, you create a random key pair (public and private key), and then send the public key to the other party. The other party uses the public key to create a shared secret and a cipher text (i.e,. the shared secret in an encrypted form). The other party sends the cipher text back to the first party. The first party uses their private key to decrypt the cipher text and get the shared secret. Now both parties have the shared secret -- this can be used as a symmetric key for encryption/decryption (e.g., with AES).
let selectedParamSet: keyof typeof Params = 'Kyber512' as const;

for (const [paramSet, path] of [
  ['Kyber512', './og-kyber/crystals-kyber-ts/src/services/kyber512.service'],
  ['Kyber768', './og-kyber/crystals-kyber-ts/src/services/kyber768.service'],
  ['Kyber1024', './og-kyber/crystals-kyber-ts/src/services/kyber1024.service'],
] as const) {
  import(path).then((imported) => {
    const KyberService = imported[`${paramSet}Service`];
    // Test on yourself
    console.log('\nTesting with', paramSet);
    console.log('\n----------------------\n');
    selectedParamSet = paramSet as keyof typeof Params;

    const Instance = new KyberService();
    const { publicKey: publicKeyU, secretKey: secretKeyU } =
      kyberCCAKEMKeyGen();
    const publicKey = Array.from(publicKeyU);
    const secretKey = Array.from(secretKeyU);

    console.log('Test: Can decrypt its own message');
    const { cipherText, sharedSecret } = kyberCCAKEMEncrypt(publicKeyU);
    const decrypted = kyberCCAKEMDecrypt(secretKeyU, cipherText);
    console.assert(
      sharedSecret.toString() === decrypted.toString(),
      'Houston, we have a bug'
    );
    console.log('Success!');

    // Sanity check that the other system works
    console.log('\nTest: the other system can decrypt its own message');
    const [c, s] = Instance.encrypt(publicKey);
    const d = Instance.decrypt(c, secretKey);
    console.assert(
      s.toString() === d.toString(),
      'the other system is not working'
    );
    console.log('Success!');

    // Other system can use our public key to encrypt a message, end we can decrypt it
    console.log(
      '\nTest: Other system can encrypt a message and we can decrypt it'
    );
    const [cipherText2, sharedSecret2] = Instance.encrypt(publicKey);
    const decrypted2 = kyberCCAKEMDecrypt(
      secretKeyU,
      new Uint8Array(cipherText2)
    );
    console.assert(
      sharedSecret2.toString() === Array.from(decrypted2).toString(),
      'could not decrypt the message correctly'
    );
    console.log('Success!');

    // We can use other system's public key to encrypt a message, and they can decrypt it
    console.log(
      '\nTest: We can encrypt a message and other system can decrypt it'
    );
    const { cipherText: cipherText3, sharedSecret: sharedSecret3 } =
      kyberCCAKEMEncrypt(new Uint8Array(publicKey));
    const decrypted3 = Instance.decrypt(Array.from(cipherText3), secretKey);
    console.assert(
      Array.from(sharedSecret3).toString() === decrypted3.toString(),
      'could not do it'
    );
    console.log('Success!');
  });
}
