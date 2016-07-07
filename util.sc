+ ArrayedCollection {

	nrOfArgs {
		var nargs = 0;
		this do: { |ia| ia.isNil.not.if { nargs = nargs + 1 } };
		^nargs;
	}

	checkNumArgs { arg lo, hi;
		var nargs = this.nrOfArgs;
		(nargs < lo).if { Error("Not enough input arguments.").throw };
		(nargs > hi).if { Error("Too many input arguments.").throw };
		^nargs
	}

	polyval { arg x;
		var n = this.size - 1, y;
		x.isKindOf(Collection).if {
			y = 0!x.size;
			x do: { |xi, i| this do: { |p, j| y[i] = y[i] + (p*(xi**(n - j))) } }
		} {
			y = 0;
			this do: { |p, i| y = y + (p*(x**(n - i))) };
		};
		^y
	}

	datawrap { arg n;
		var nx, nr, a;

		(n < 1).if { Error("n needs to be a positive integer larger than 0").throw };
		(n <= 1).if { ^this.sum };
		(this.size == n).if { ^this };

		nx = this.size;

		nr = (this.size/n).ceil.asInteger;
		a = this.extend((nr*n), 0);   // zero pad
		(n > nx).if {
			^a
		} {
			^a.reshape(nr, n).flop.performUnaryOp(\sum)
		}
	}

}

+ Symbol {

	toLower { ^this.asString.toLower.asSymbol }

}