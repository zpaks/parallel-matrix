/* this is a distributed program run in parallel over multiple locales.
 */

use Random, Time;
use BlockDist, CyclicDist;

/*
  Here we set up several config consts that represent, in order:
    - ``N``, the dimension for the square array ``A``
    - ``K``, the second dimension for arrays ``X`` and ``B``
    - ``seed``, a seed for random number generation
*/

config const N = 400;
config const K = 400;

config const SQRTHREADS = 2;

config const seed = 100;
config const seed2 = 100;

// create a variable to create the space needed in block distribution
const Space = {1..N, 1..K};

// Create the arrays ``A``, ``X``, and ``B``. Fill ``A`` and ``X`` with random
// values.

var A = blockDist.createArray({1..N, 1..N}, real);
fillRandom(A, seed);

var X = blockDist.createArray({1..N, 1..N}, real);
fillRandom(X, seed2);

var B = blockDist.createArray({1..N, 1..N}, real);

var AA = blockDist.createArray({1..N*N}, real);
fillRandom(AA, seed);

var XX = blockDist.createArray({1..N*N}, real);
fillRandom(XX, seed2);

var BB = blockDist.createArray({1..N*N}, real);


// variables for the time taken by the program.
var t:Timer;


// Matrix multiply ``A*X``, store result in ``B``
proc Mmultiply(N: int, K: int): void
{
    t.start(); // start timer
    forall i in 1..N do {
      forall j in 1..K do { 
        forall k in 1..N do { 
           B[i,j] += A[i,k] * X[k,j];
        }
      }
    }
    t.stop(); // stop timer
    writeln(t.elapsed()," seconds elapsed");
    t.clear();
}

// Matrix multiply ``A*X``, store result in ``B``
proc MmultiplySqrt(N: int, K: int): void
{
    t.start(); // start timer
    
    var blockSize: int = N / SQRTHREADS;
    var threads = SQRTHREADS * SQRTHREADS;
    forall p in 1..threads do {
        var rowIndex: int = (p - 1) / SQRTHREADS;
        var columnIndex: int = (p - 1) % SQRTHREADS;

        for ri in 1..blockSize do {
            for rj in 1..blockSize do {
                var i: int = rowIndex * blockSize + ri;
                var j: int = columnIndex * blockSize + rj;

                var sum: real = 0;
                for k in 1..N do { 
                    sum += A[i,k] * X[k,j];
                    // sum += 1232;
                }

                B[i,j] = sum;
            }
        }
    }

    t.stop(); // stop timer
    writeln(t.elapsed()," seconds elapsed");
    t.clear();
}

// Matrix multiply ``A*X``, store result in ``B``
proc MmultiplySqrtV2(N: int, K: int): void
{
    t.start(); // start timer
    
    var block_num = SQRTHREADS;
    var m_size = N;
    block_num = if block_num > m_size then m_size else block_num;
    var block_size = m_size / block_num;

    for offset in 0..(block_num - 1) {
        forall i in 0..(block_num * block_num - 1) {
            var j = (i / block_num + i + offset) % block_num;

            var la = if (i / block_num == block_num - 1) then m_size - block_size * (block_num - 1) else block_size;
            var lb = if (i % block_num == block_num - 1) then m_size - block_size * (block_num - 1) else block_size;
            var lc = if (j == block_num - 1) then m_size - block_size * (block_num - 1) else block_size;

            var ma = block_size * (i / block_num) * m_size + j * block_size;
            var mb = block_size * (i % block_num) + j * block_size * m_size;
            var mc = block_size * (i / block_num) * m_size + block_size * (i % block_num);

            // writeln(ma, mb, mc);

            for ii in 0..la-1 do {
                for jj in 0..lb-1 do  {
                    var wyn :real = 0;
                    for kk in 0..lc-1 do {
                        // wyn += AA(ii, kk) * XX(kk, jj);
                        wyn += AA[ma + ii * m_size + kk + 1] * XX[mb + kk * m_size + jj + 1];
                    }
                    // C(i, j) += wyn;
                    BB[mc + ii * m_size + jj + 1] += wyn;
                }
            }

			// multiplyMatrix(a + ia, la, b + ib, lb, temp + mc, lc, m_size);
		}
    }

    t.stop(); // stop timer
    writeln(t.elapsed()," seconds elapsed");
    t.clear();
}

// loop to print results of 10 executions of Mmultiply proc. 
for f in 1..5 do {
    MmultiplySqrtV2(N, K);
 }

writeln(" 10 calculations complete");
