import { SHA3, SHAKE } from 'sha3';

// TODO: Define param sets
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

// This is the rejection sampling
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

      // If we didn't get enough bytes, we need to sample more
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

const QINV = -3327 as const; // q ** -1 & (2**16)
function montgomeryReduce(a: number): number {
  const t = a * QINV;
  return (a = t * Params[selectedParamSet].q) >> 16;
}

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
          const t = montgomeryReduce(zeta * r[j + len]) >>> 0;
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
    result[4 * i] = montgomeryReduce(a[4 * i + 1] * b[4 * i + 1]);
    result[4 * i] = montgomeryReduce(result[4 * i] + zeta);
    result[4 * i] += montgomeryReduce(a[4 * i] * b[4 * i]);
    result[4 * i + 1] = montgomeryReduce(a[4 * i] + b[4 * i + 1]);
    result[4 * i + 1] += montgomeryReduce(a[4 * i + 1] + b[4 * i]);

    result[4 * i + 2] = montgomeryReduce(a[4 * i + 3] * b[4 * i + 3]);
    result[4 * i + 2] = montgomeryReduce(result[4 * i + 2] + -zeta);
    result[4 * i + 2] += montgomeryReduce(a[4 * i + 2] * b[4 * i + 2]);
    result[4 * i + 3] = montgomeryReduce(a[4 * i + 2] + b[4 * i + 3]);
    result[4 * i + 3] += montgomeryReduce(a[4 * i + 3] + b[4 * i + 2]);
  }
  return result;
}

// Add two polynomials. Mutates target. Assume target, a, and b are the same length
function addPolynomials(
  target: Uint16Array,
  a: Uint16Array,
  b: Uint16Array
): void {
  for (let i = 0; i < Params[selectedParamSet].n; i++) {
    target[i] = a[i] + b[i];
  }
}

// Add two polynomials. Mutates target. Assume target, a, and b are the same length
function addPolynomialVectors(
  target: Uint16Array[],
  a: Uint16Array[],
  b: Uint16Array[]
): void {
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    addPolynomials(target[i], a[i], b[i]);
  }
}

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

function applyBarrettReductionToVectorInPlace(r: Uint16Array[]): void {
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    applyBarrettReductionInPlace(r[i]);
  }
}

function applyMatrix(
  A: Uint16Array[][],
  s: Uint16Array[],
  resultLen: number
): Uint16Array[] {
  const result: Uint16Array[] = new Array(resultLen);
  const f = (1 << 32) % Params[selectedParamSet].q;
  for (let i = 0; i < Params[selectedParamSet].k; i++) {
    result[i] = multiplyPolynomials(A[i][0], s[0]);
    for (let j = 1; j < Params[selectedParamSet].k; j++) {
      const t = multiplyPolynomials(A[i][j], s[j]);
      addPolynomials(result[i], result[i], t);
    }

    applyBarrettReductionInPlace(result[i]);
    for (let j = 0; j < Params[selectedParamSet].n; j++) {
      result[i][j] = montgomeryReduce(result[i][j] * f);
    }
  }
  return result;
}

// Output of length 384 bytes
const POLYBYTES = 384 as const;
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

// Generate key pair. Output is a tuple of secret key and public key.
const SYMBYTES = 32 as const;
export function kyberCPAPKEKeyGen() {
  // Generate random bytes for seed
  const randomBytes = Buffer.alloc(SYMBYTES);
  for (let i = 0; i < SYMBYTES; i++) {
    randomBytes[i] = Math.floor(Math.random() * 256); // TODO: Requires a safer random number generator
  }
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

  const t = applyMatrix(A, s, Params[selectedParamSet].k);
  addPolynomialVectors(t, t, e);
  applyBarrettReductionToVectorInPlace(t);

  const polyVecBytes = Params[selectedParamSet].k * POLYBYTES;
  const publicKey = new Uint8Array(polyVecBytes + SYMBYTES);
  publicKey.set(encode(t));
  publicKey.set(publicSeed, polyVecBytes);
  const secretKey = encode(s);

  return { secretKey, publicKey };
}

// TODO: Exchange keys

// TODO: Encrypt message with public key

// TODO: Decrypt message with private key

// TODO: Implement CLI
const selectedParamSet: keyof typeof Params = 'Kyber512' as const;

console.log(kyberCPAPKEKeyGen());
