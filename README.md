# CZT
Chirp z-transform implementation for SuperCollider

## Example

// Make a test signal
l = 73;
x = l collect: { |n| 1 + cos(8*n*2pi/l) + cos(3*n*2pi/l) };

// Do the transform (equivalent to an arbitrary size FFT, no need for power of two FFT size!!!)
y = x.czt;
y.abs.plot(discrete: true);
