// Ported to golang from https://getkerf.wordpress.com/2016/03/30/the-best-algorithm-no-one-knows-about/
package vitter

import (
	"math"
	"math/rand"
	"time"
)

func round(val float64) int {
	if val < 0 {
		return int(val - 0.5)
	}
	return int(val + 0.5)
}

func random_double() float64 {
	//[0,1) uniformly random double
	r := rand.New(rand.NewSource(time.Now().UnixNano()))
	return r.Float64()
}

//Vitter, J.S. - An Efficient Algorithm for Sequential Random Sampling - ACM Trans. Math. Software 11 (1985), 37-57.
//'a' is space allocated for the hand
//'n' is the size of the hand
//'N' is the upper bound on the random card values
func Vitter(a []int64, n, N int64) []int64 {

	var i, j, t, qu1, S, negalphainv, threshold int64

	j = -1
	qu1 = -n + 1 + N
	negalphainv = -13
	threshold = -negalphainv * n

	var nreal, Nreal, ninv, nmin1inv, Vprime float64

	nreal = float64(n)
	Nreal = float64(N)
	ninv = float64(1.0 / n)
	nmin1inv = float64(1.0 / (n - 1))
	Vprime = math.Exp(math.Log(random_double()) * ninv)

	var qu1real, negSreal, U, X, y1, y2, top, bottom, limit float64

	qu1real = -nreal + 1.0 + Nreal

	for n > 1 && threshold < N {
		nmin1inv = 1.0 / (-1.0 + nreal)

		for {
			for {
				X = Nreal * (-Vprime + 1.0)
				S = int64(math.Floor(X))

				if S < qu1 {
					break
				}

				Vprime = math.Exp(math.Log(random_double()) * ninv)
			}

			U = random_double()

			negSreal = float64(-S)

			y1 = math.Exp(math.Log(U*Nreal/qu1real) * nmin1inv)

			Vprime = y1 * (-X/Nreal + 1.0) * (qu1real / (negSreal + qu1real))

			if Vprime <= 1.0 {
				break
			}

			y2 = 1.0

			top = -1.0 + Nreal

			if -1+n > S {
				bottom = -nreal + Nreal
				limit = float64(-S + N)
			} else {
				bottom = -1.0 + negSreal + Nreal
				limit = float64(qu1)
			}

			t = N - 1
			for float64(t) >= limit {
				y2 = (y2 * top) / bottom
				top--
				bottom--
				t--
			}

			if Nreal/(-X+Nreal) >= y1*math.Exp(math.Log(y2)*nmin1inv) {
				Vprime = math.Exp(math.Log(random_double()) * nmin1inv)
				break
			}

			Vprime = math.Exp(math.Log(random_double()) * ninv)
		}

		j += S + 1

		a[i] = j
		i++

		N = -S + (-1 + N)

		Nreal = negSreal + (-1.0 + Nreal)

		n--
		nreal--
		ninv = nmin1inv

		qu1 = -S + qu1
		qu1real = negSreal + qu1real

		threshold += negalphainv
	}

	if n > 1 {
		vitter_a(a[i:], n, N, j) // if i>0 then n has been decremented
	} else {
		S = int64(math.Floor(float64(N) * Vprime))

		j += S + 1

		a[i] = j
		i++
	}

	return a
}

func vitter_a(a []int64, n, N, j int64) []int64 {
	var S, i int64
	var V, quot float64
	top := float64(N - n)
	Nreal := float64(N)


	for n >= 2 {
		V = random_double()
		S = 0
		quot = top / Nreal

		for quot > V {
			S++
			top--
			Nreal--
			quot = (quot * top) / Nreal
		}

		j += S + 1

		a[i] = j
		i++

		Nreal--

		n--
	}

	S = int64(math.Floor(float64(round(Nreal)) * random_double()))

	j += S + 1

	a[i] = j
	i++

	return a
}
