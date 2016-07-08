+ ArrayedCollection {

	/**
	 * Digital filter frequency response.
     * [h, w] = freqz(b, a, n) returns the n-point complex frequency response
	 * vector h and the n-point frequency vector w in radians/sample of
	 * the filter:
     *           jw               -jw              -jmw
	 *    jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
	 * H(e) = ---- = ------------------------------------
	 *           jw               -jw              -jnw
	 *        A(e)    a(1) + a(2)e + .... + a(n+1)e
	 * given numerator and denominator coefficients in vectors b and a. The
	 * frequency response is evaluated at N points equally spaced around the
	 * upper half of the unit circle. If N isn't specified, it defaults to
     * 512.
	 *
	 * NOTE: this is a SC port of MATLAB's freqz function
	 */
	freqz { arg n = 512, fs = nil, w = nil, range = \onesided, plot = true;
		var nargin, fvflag = false, centerdc = false, b, a, h;

		(this.shape.size == 2 and: { this.shape[0] == 2 }).if {
			b = this[0].copy;
			a = this[1].copy;
		} {
			(this.shape.size == 1).if {
				b = this.copy;
				a = 1;
			} {
				Error("freqz expects either a single array (FIR) or an array of two arrays (IIR)").throw
			}
		};

		#[\onesided, \twosided].includes(range).not.if {
			Error("Range should either be 'onesided' or 'twosided'.").throw
		};

		w.notNil.if {
			w.isKindOf(Collection).not.if { Error("w/f needs to be an array of frequency points.").throw };
			(w.size < 2).if { Error("w/f must contain at least two frequency points.").throw };
			fvflag = true;
			n.notNil.if { "The size of the frequency vector w/f is used to determine the number of points in the complex frequency response.".warn };
		};

		(a.size == 0).if {
			#h, w = this.prFIRFreqz(b/a, w, fs, n, fvflag, range);
		} {
			#h, w = this.prIIRFreqz(b, a, w, fs, n, fvflag, range);
		};

		plot.if {
			[h, w].prFreqzPlot();
		}

		^[h, w]
	}

	prFIRFreqz { arg b, w, fs, nfft, fvflag, range;
		var n, digw, s, h;

		n = b.size;

		fvflag.if {
			/**
			 * Frequency vector specified. Use Horner's method of polynomial
			 * evaluation at the frequency points and divide the numerator
			 * by the denominator.
			 */
			fs.notNil.if {
				digw = 2pi*w/fs;   // convert from Hz to rad/sample for computational purposes
			} {
				digw = w;
			};

			s = exp(Complex(0, 1)*digw);
			h = b.polyval(s)/exp(Complex(0, 1)*digw*(n - 1));
		} {
			// freqvector not specified, use nfft and range in calculation
			s = #[\twosided, \onesided].indexOf(range) + 1;

			(s*nfft < n).if {
				// data is larger than FFT points, wrap modulo s*nfft
				b = b.datawrap(s*nfft);
			};

			h = b.czt(s*nfft);
			// when range = 'onesided', we computed a 2*nfft point FFT, now we talk half the result
			h = h[0..nfft - 1];
			w = this.prFreqzFreqVec(nfft, fs, s);
		};

		^[h, w]
	}

	prIIRFreqz { arg b, a, w, fs, nfft, fvflag, range;
		var nb, na, n, digw, s, h;

		nb = b.size;
		na = a.size;
		(nb > na).if { a = a ++ (0!(nb - na)) };  // make a and b the same length
		(na > nb).if { b = b ++ (0!(na - nb)) };
		n = a.size;   // this will be the new length of both num and den

		fvflag.if {
			/**
			 * Frequency vector specified. Use Horner's method of polynomial
			 * evaluation at the frequency points and divide the numerator
			 * by the denominator.
			 */
			fs.notNil.if {
				digw = 2pi*w/fs;   // convert from Hz to rad/sample for computational purposes
			} {
				digw = w;
			};

			s = exp(Complex(0, 1)*digw);
			h = b.polyval(s)/a.polyval(s);
		} {
			// freqvector not specified, use nfft and range in calculation
			s = #[\twosided, \onesided].indexOf(range) + 1;

			(s*nfft < n).if {
				// data is larger than FFT points, wrap modulo s*nfft
				b = b.datawrap(s*nfft);
				a = a.datawrap(s*nfft);
			};

			h = b.czt(s*nfft)/a.czt(s*nfft);
			// when range = 'onesided', we computed a 2*nfft point FFT, now we talk half the result
			h = h[0..nfft - 1];
			w = this.prFreqzFreqVec(nfft, fs, s);
		};

		^[h, w]
	}

	prFreqzFreqVec { arg nfft, fs, s;
		var w, deltaF;

		fs = fs ? 2pi;
		s = s ? 2;

		switch(s,
			1, {   // 0-2pi
				deltaF = fs/nfft;
				w = (0..nfft - 1)/(nfft - 1)*(fs - deltaF);

				// There can still be some minor round off errors. Fix the known points,
				// those near pi and 2pi.
				(nfft%2).booleanValue.if {
					w[(nfft + 1).div(2) - 1] = 0.5*fs*(1 - nfft.reciprocal);
					w[(nfft + 1).div] = 0.5*fs*(1 + nfft.reciprocal);
				} {
					// make sure we hit fs/2 exactly for the 1/2 Nyquist point
					w[nfft.div(2)] = 0.5*fs;
				};
				w[nfft - 1] = fs - (fs/nfft);
			},
			2, {   // 0-pi
				deltaF = 0.5*fs/nfft;
				w = (0..nfft - 1)/(nfft - 1)*(0.5*fs - deltaF);

				w[nfft - 1] = 0.5*fs*(1 - nfft.reciprocal);
			},
			3, {   // -pi-pi
				var wmin, wmax;
				deltaF = fs/nfft;

				(nfft%2).booleanValue.if {   // ODD, don't include Nyquist
					wmin = -0.5*(fs - deltaF);
					wmax = 0.5*(fs - deltaF);
				} {   // EVEN, include Nyquist point in the negative freq
					wmin = -0.5*fs;
					wmax = 0.5*fs - deltaF;
				};
				w = wmin + ((0..nfft - 1)/(nfft - 1)*(wmax - wmin));
				(nfft%2).booleanValue.if {   // ODD
					w[(nfft + 1).div(2) - 1] = 0;
				} {
					w[nfft.div(2)] = 0;
				}
			}
		);

		^w
	}

	prFreqzPlot {
	}

}