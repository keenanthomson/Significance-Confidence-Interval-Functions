function cdfNormal (x, mean, standardDeviation) {
  return (1 - mathjs.erf((mean - x ) / (Math.sqrt(2) * standardDeviation))) / 2
}

function erf (x) {
  const y = Math.abs(x)

  if (y >= MAX_NUM) {
    return Math.sign(x)
  }
  if (y <= THRESH) {
    return Math.sign(x) * erf1(y)
  }
  if (y <= 4.0) {
    return Math.sign(x) * (1 - erfc2(y))
  }
  return Math.sign(x) * (1 - erfc3(y))
}

function erf1 (y) {
  const ysq = y * y
  let xnum = P[0][4] * ysq
  let xden = ysq
  let i

  for (i = 0; i < 3; i += 1) {
    xnum = (xnum + P[0][i]) * ysq
    xden = (xden + Q[0][i]) * ysq
  }
  return y * (xnum + P[0][3]) / (xden + Q[0][3])
}

function erfc2 (y) {
  let xnum = P[1][8] * y
  let xden = y
  let i

  for (i = 0; i < 7; i += 1) {
    xnum = (xnum + P[1][i]) * y
    xden = (xden + Q[1][i]) * y
  }
  const result = (xnum + P[1][7]) / (xden + Q[1][7])
  const ysq = parseInt(y * 16) / 16
  const del = (y - ysq) * (y + ysq)
  return Math.exp(-ysq * ysq) * Math.exp(-del) * result
}

function erfc3 (y) {
  let ysq = 1 / (y * y)
  let xnum = P[2][5] * ysq
  let xden = ysq
  let i

  for (i = 0; i < 4; i += 1) {
    xnum = (xnum + P[2][i]) * ysq
    xden = (xden + Q[2][i]) * ysq
  }
  let result = ysq * (xnum + P[2][4]) / (xden + Q[2][4])
  result = (SQRPI - result) / y
  ysq = parseInt(y * 16) / 16
  const del = (y - ysq) * (y + ysq)
  return Math.exp(-ysq * ysq) * Math.exp(-del) * result
}

const THRESH = 0.46875

const SQRPI = 5.6418958354775628695e-1

const MAX_NUM = Math.pow(2, 53);

const P = [[
  3.16112374387056560e00, 1.13864154151050156e02,
  3.77485237685302021e02, 3.20937758913846947e03,
  1.85777706184603153e-1
], [
  5.64188496988670089e-1, 8.88314979438837594e00,
  6.61191906371416295e01, 2.98635138197400131e02,
  8.81952221241769090e02, 1.71204761263407058e03,
  2.05107837782607147e03, 1.23033935479799725e03,
  2.15311535474403846e-8
], [
  3.05326634961232344e-1, 3.60344899949804439e-1,
  1.25781726111229246e-1, 1.60837851487422766e-2,
  6.58749161529837803e-4, 1.63153871373020978e-2
]]

const Q = [[
  2.36012909523441209e01, 2.44024637934444173e02,
  1.28261652607737228e03, 2.84423683343917062e03
], [
  1.57449261107098347e01, 1.17693950891312499e02,
  5.37181101862009858e02, 1.62138957456669019e03,
  3.29079923573345963e03, 4.36261909014324716e03,
  3.43936767414372164e03, 1.23033935480374942e03
], [
  2.56852019228982242e00, 1.87295284992346047e00,
  5.27905102951428412e-1, 6.05183413124413191e-2,
  2.33520497626869185e-3
]]

// console.log(erf(0.2));

/*
inverse cumulative distribution function/helpers (quantile function)
*/

  function erfcinv(p) {
    var j = 0;
    var x, err, t, pp;
    if (p >= 2) {
      return -100;
    }
    if (p <= 0) {
      return 100;
    }
    pp = (p < 1) ? p : 2 - p;
    t = Math.sqrt(-2 * Math.log(pp / 2));
    x = -0.70711 * ((2.30753 + t * 0.27061) /
      (1 + t * (0.99229 + t * 0.04481)) - t);
    for (; j < 2; j++) {
      err = erfc(x) - pp;
      x += err / (1.12837916709551257 * Math.exp(-x * x) - x * err);
    }
    return (p < 1) ? x : -x;
  }

  function erfc(x) {
    return 1 - erf(x);
  }

  function erf(x) {
    var cof = [-1.3026537197817094, 6.4196979235649026e-1, 1.9476473204185836e-2,
      -9.561514786808631e-3, -9.46595344482036e-4, 3.66839497852761e-4,
      4.2523324806907e-5, -2.0278578112534e-5, -1.624290004647e-6,
      1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8,
      6.529054439e-9, 5.059343495e-9, -9.91364156e-10,
      -2.27365122e-10, 9.6467911e-11, 2.394038e-12,
      -6.886027e-12, 8.94487e-13, 3.13092e-13,
      -1.12708e-13, 3.81e-16, 7.106e-15,
      -1.523e-15, -9.4e-17, 1.21e-16,
      -2.8e-17
    ];
    var j = cof.length - 1;
    var isneg = false;
    var d = 0;
    var dd = 0;
    var t, ty, tmp, res;

    if (x < 0) {
      x = -x;
      isneg = true;
    }

    t = 2 / (2 + x);
    ty = 4 * t - 2;

    for (; j > 0; j--) {
      tmp = d;
      d = ty * d - dd + cof[j];
      dd = tmp;
    }

    res = t * Math.exp(-x * x + 0.5 * (cof[0] + ty * d) - dd);
    return isneg ? res - 1 : 1 - res;
  }

  function inv(p, mean, std) {
    return -1.41421356237309505 * std * erfcinv(2 * p) + mean;
  }

