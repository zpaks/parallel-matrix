use Random, Time;
use BlockDist, CyclicDist;

config const n = 10000000;
config const threads = 1;

var seed: int = 589494289 ; // seed for random number generator
writeln("Number of points = ", n);
writeln("Random number seed = ", seed);
var A = blockDist.createArray({1..threads}, int);

var t: Timer;
t.start(); // start timer

forall i in 1..threads do {
    var rs = new RandomStream(real, seed, parSafe=false);
    var count = 0;
    var part = n / threads;
    for i in 1..part do {
        if (rs.getNext()**2 + rs.getNext()**2) <= 1.0 then
            count += 1;
    }
    A[i] = count;
}

var cnt = 0;
for i in 1..threads do {
    cnt += A[i];
}

t.stop(); // stop timer

writeln(cnt * 4.0 / n);
writeln(t.elapsed()," seconds elapsed");
