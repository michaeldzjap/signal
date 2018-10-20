# signal
Port of some MATLAB functions from the signal toolbox to SuperCollider

## Examples

### Chirp z-transform
```
// Make a test signal
l = 73;
x = l collect: { |n| 1 + cos(8*n*2pi/l) + cos(3*n*2pi/l) };

// Do the transform (equivalent to an arbitrary size FFT, no need for power of two FFT size!!!)
y = x.czt;
y.abs.plot("Magnitude of y(n)", Rect(0, 0, 800, 400), true);
```

### Parks/McClellan optimal FIR filter design
```
// lowpass example
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.05, 0.1], [1, 0], [0.01, 0.001]);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n, f, a, w);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot("Frequency response H(z)", Rect(0, 0, 800, 400), numChannels: 2);
)

// bandpass example
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.2, 0.25, 0.6, 0.7], [0, 1, 0], [0.001, 0.01, 0.01]);

// fix bump in magnitude response (comment out this line to see the difference)
#f, a, w = FIR.convertToConstraint(f, [0, 1, 0], [0.001, 0.01, 0.01], 0.0005, \type_a);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n, f, a, w);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot("Frequency response H(z)", Rect(0, 0, 800, 400), numChannels: 2);
)
```
