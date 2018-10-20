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

// highpass example
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.02, 0.05], [0, 1], [0.001, 0.01]);

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

// bandstop example
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.15, 0.3, 0.6, 0.65], [1, 0, 1], [0.01, 0.001, 0.01]);

// fix bump in magnitude response (comment out this line to see the difference)
#f, a, w = FIR.convertToConstraint(f, [1, 0, 1], [0.01, 0.001, 0.01], 0.0005, \type_a);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n, f, a, w);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot("Frequency response H(z)", Rect(0, 0, 800, 400), numChannels: 2);
)

// five-band filter
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.17, 0.23, 0.47, 0.53, 0.67, 0.73, 0.82, 0.88], [0, 1, 0, 1, 0], [0.001, 0.01, 0.001, 0.01, 0.001]);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n + 4, f, a, w, impType: 1);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot("Frequency response H(z)", Rect(0, 0, 800, 400), numChannels: 2);
)

// another five-band filter
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.1, 0.15, 0.3, 0.35, 0.75, 0.8, 0.85, 0.9], [1, 0, 1, 0, 1], [0.01, 0.001, 0.01, 0.001, 0.01]);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n + 4, f, a, w, impType: 1);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot("Frequency response H(z)", Rect(0, 0, 800, 400), numChannels: 2);
)

/**
 * Example of an L-th band (Nyquist) lowpass FIR filter that can be used for 4x
 * oversampling (L = 4 in this case). Note how every 4-th coefficient is 0 so
 * that every 4-th tap of a polyphase implementation basically reduces from a
 * convolution operation to a simple scaling of the relevant input sample.
 * The frequency range [0 - 1] for FIR.parksMcClellan corresponds to the range
 * [0 - fs / 2] Hz (or equivalently, [0 - pi] radians).
 *
 * The filter length for an L-th band FIR has to be uneven (hence n needs to be
 * even).
 */
(
var n = 510, f = [0, 0.24, 0.26, 1], a = [1, 1, 0, 0], w = [1, 10], h, hf, hr, err, p;

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n, f, a, w);
#hf, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [hf.abs.max(1e-06).ampdb, hf.phase.unwrap].lace(hf.size * 2).plot("Frequency response H(z)", Rect(0, 0, 800, 400), numChannels: 2);
p.plots[0].domainSpec = [0, 0.5].asSpec; // 0 - fs / 2
p.plots[1].domainSpec = [0, 1].asSpec;
p.refresh;

q = h[n.div(2.5)..h.lastIndex - n.div(2.5)].plot("Impulse response h(n)", Rect(0, 500, 800, 300), true);
q.domainSpecs = [0, n - n.div(2.5) - n.div(2.5)].asSpec;
q.refresh;
)
