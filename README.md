# LightECC

<div align="center">

[![PyPI Downloads](https://static.pepy.tech/personalized-badge/lightecc?period=total&units=international_system&left_color=grey&right_color=blue&left_text=downloads)](https://pepy.tech/project/lightecc)
[![Stars](https://img.shields.io/github/stars/serengil/LightECC?color=yellow&style=flat&label=%E2%AD%90%20stars)](https://github.com/serengil/LightECC/stargazers)
[![Tests](https://github.com/serengil/LightECC/actions/workflows/tests.yml/badge.svg)](https://github.com/serengil/LightECC/actions/workflows/tests.yml)
[![License](http://img.shields.io/:license-MIT-green.svg?style=flat)](https://github.com/serengil/LightECC/blob/master/LICENSE)

[![Blog](https://img.shields.io/:blog-sefiks.com-blue.svg?style=flat&logo=wordpress)](https://sefiks.com)
[![YouTube](https://img.shields.io/:youtube-@sefiks-red.svg?style=flat&logo=youtube)](https://www.youtube.com/@sefiks?sub_confirmation=1)
[![Twitter](https://img.shields.io/:follow-@serengil-blue.svg?style=flat&logo=x)](https://twitter.com/intent/user?screen_name=serengil)

[![Support me on Patreon](https://img.shields.io/endpoint.svg?url=https%3A%2F%2Fshieldsio-patreon.vercel.app%2Fapi%3Fusername%3Dserengil%26type%3Dpatrons&style=flat)](https://www.patreon.com/serengil?repo=lightecc)
[![GitHub Sponsors](https://img.shields.io/github/sponsors/serengil?logo=GitHub&color=lightgray)](https://github.com/sponsors/serengil)
[![Buy Me a Coffee](https://img.shields.io/badge/-buy_me_a%C2%A0coffee-gray?logo=buy-me-a-coffee)](https://buymeacoffee.com/serengil)

</div>

<p align="center"><img src="https://raw.githubusercontent.com/serengil/LightECC/master/images/starfish.jpg" width="240" height="240"></p>

LightECC is a lightweight elliptic curve cryptography library for its arithmetic for python. It is a hybrid library wrapping many elliptic curve forms such as [Weierstrass](https://sefiks.com/2016/03/13/the-math-behind-elliptic-curve-cryptography/), [Koblitz](https://sefiks.com/2016/03/13/the-math-behind-elliptic-curves-over-binary-field/) and [Edwards](https://sefiks.com/2018/12/19/a-gentle-introduction-to-edwards-curves/).

# Elliptic Curve Arithmetic

Building an elliptic curve cryptosystem is very straightforward in LightECC. You basically need to initialize a LightECC object with a form name and a curve name. By default, it constructs elliptic curves in Weierstras form.After that, you can retrieve the base point of the curve and perform various elliptic curve arithmetic operations, including addition, subtraction, multiplication, and division.

```python
from lightecc import LightECC

forms = ["weierstrass", "koblitz", "edwards"]

ec = LightECC(
    form_name = "edwards",
    curve_name = "ed25519",
)

# get the base point
G = ec.G

# addition
_2G = G + G
_3G = _2G + G
_5G = _3G + _2G
_10G = _5G + _5G

# subtraction
_9G = _10G - G

# multiplication
_20G = 20 * G
_50G = 50 * G

# division
_25G = _50G / G
```

Here, the [double-and-add](https://sefiks.com/2016/03/27/double-and-add-method/) method is used for multiplication, allowing it to be performed very quickly, regardless of the size of the multiplier.

On the other hand, division requires solving the [elliptic curve discrete logarithm problem](https://sefiks.com/2018/02/28/attacking-elliptic-curve-discrete-logarithm-problem/), which is computationally difficult.

# Point at Infinity or Neutral & Identity Element

The order of the elliptic curve is defined by the argument n in the constructed LightECC object. This represents the total number of points on the curve. It also serves as the [neutral or identity element](https://sefiks.com/2023/09/29/understanding-identity-element-in-elliptic-curves/) of the curve, meaning that adding this point to any other point does not change the result. Additionally, elliptic curves exhibit cyclic group properties beyond this point.

```python
ec = LightECC()

# order of elliptic curve
n = ec.n

# neutral element
nG = n * G

# scalar multiplication
_17G = 17 * G

# proof of work for neutralism
assert _17G == _17G + nG

# proof of work for cyclic group
assert (n + 1) * G == G
assert (n + 2) * G == 2 * G
```

# Supported Curves

Below is a list of elliptic curves supported by LightECC. Each curve has a specific order (n), which defines the number of points in the finite field. The order directly impacts the cryptosystem's security strength. A higher order typically corresponds to a stronger cryptosystem, making it more resistant to cryptographic attacks.

## Edwards Curves

| form | curve | field | n (bits) |
| --- | --- | --- | --- |
| edwards | e521 | prime | 519 |
| edwards | id-tc26-gost-3410-2012-512-paramsetc | prime | 510 |
| edwards | numsp512t1 | prime | 510 |
| edwards | ed448 | prime | 446 |
| edwards | curve41417 | prime | 411 |
| edwards | numsp384t1 | prime | 382 |
| edwards | id-tc26-gost-3410-2012-256-paramseta | prime | 255 |
| edwards | ed25519 | prime | 254 |
| edwards | mdc201601 | prime | 254 |
| edwards | numsp256t1 | prime | 254 |
| edwards | jubjub | prime | 252 |

## Weierstass Form

| form | curve | field | n (bits) |
| --- | --- | --- | --- |
| weierstrass | bn638 | prime | 638 |
| weierstrass | bn606 | prime | 606 |
| weierstrass | bn574 | prime | 574 |
| weierstrass | bn542 | prime | 542 |
| weierstrass | p521 | prime | 521 |
| weierstrass | brainpoolp512r1 | prime | 512 |
| weierstrass | brainpoolp512t1 | prime | 512 |
| weierstrass | fp512bn | prime | 512 |
| weierstrass | numsp512d1 | prime | 512 |
| weierstrass | gost512 | prime | 511 |
| weierstrass | bn510 | prime | 510 |
| weierstrass | bn478 | prime | 478 |
| weierstrass | bn446 | prime | 446 |
| weierstrass | bls12-638 | prime | 427 |
| weierstrass | bn414 | prime | 414 |
| weierstrass | brainpoolp384r1 | prime | 384 |
| weierstrass | brainpoolp384t1 | prime | 384 |
| weierstrass | fp384bn | prime | 384 |
| weierstrass | numsp384d1 | prime | 384 |
| weierstrass | p384 | prime | 384 |
| weierstrass | bls24-477 | prime | 383 |
| weierstrass | bn382 | prime | 382 |
| weierstrass | curve67254 | prime | 380 |
| weierstrass | bn350 | prime | 350 |
| weierstrass | brainpoolp320r1 | prime | 320 |
| weierstrass | brainpoolp320t1 | prime | 320 |
| weierstrass | bn318 | prime | 318 |
| weierstrass | bls12-455 | prime | 305 |
| weierstrass | bls12-446 | prime | 299 |
| weierstrass | bn286 | prime | 286 |
| weierstrass | brainpoolp256r1 | prime | 256 |
| weierstrass | brainpoolp256t1 | prime | 256 |
| weierstrass | fp256bn | prime | 256 |
| weierstrass | gost256 | prime | 256 |
| weierstrass | numsp256d1 | prime | 256 |
| weierstrass | p256 | prime | 256 |
| weierstrass | secp256k1 | prime | 256 |
| weierstrass | tom256 | prime | 256 |
| weierstrass | bls12-381 | prime | 255 |
| weierstrass | pallas | prime | 255 |
| weierstrass | tweedledee | prime | 255 |
| weierstrass | tweedledum | prime | 255 |
| weierstrass | vesta | prime | 255 |
| weierstrass | bn254 | prime | 254 |
| weierstrass | fp254bna | prime | 254 |
| weierstrass | fp254bnb | prime | 254 |
| weierstrass | bls12-377 | prime | 253 |
| weierstrass | curve1174 | prime | 249 |
| weierstrass | mnt4 | prime | 240 |
| weierstrass | mnt5-1 | prime | 240 |
| weierstrass | mnt5-2 | prime | 240 |
| weierstrass | mnt5-3 | prime | 240 |
| weierstrass | prime239v1 | prime | 239 |
| weierstrass | prime239v2 | prime | 239 |
| weierstrass | prime239v3 | prime | 239 |
| weierstrass | secp224k1 | prime | 225 |
| weierstrass | brainpoolp224r1 | prime | 224 |
| weierstrass | brainpoolp224t1 | prime | 224 |
| weierstrass | curve4417 | prime | 224 |
| weierstrass | fp224bn | prime | 224 |
| weierstrass | p224 | prime | 224 |
| weierstrass | bn222 | prime | 222 |
| weierstrass | curve22103 | prime | 218 |
| weierstrass | brainpoolp192r1 | prime | 192 |
| weierstrass | brainpoolp192t1 | prime | 192 |
| weierstrass | p192 | prime | 192 |
| weierstrass | prime192v2 | prime | 192 |
| weierstrass | prime192v3 | prime | 192 |
| weierstrass | secp192k1 | prime | 192 |
| weierstrass | bn190 | prime | 190 |
| weierstrass | secp160k1 | prime | 161 |
| weierstrass | secp160r1 | prime | 161 |
| weierstrass | secp160r2 | prime | 161 |
| weierstrass | brainpoolp160r1 | prime | 160 |
| weierstrass | brainpoolp160t1 | prime | 160 |
| weierstrass | mnt3-1 | prime | 160 |
| weierstrass | mnt3-2 | prime | 160 |
| weierstrass | mnt3-3 | prime | 160 |
| weierstrass | mnt2-1 | prime | 159 |
| weierstrass | mnt2-2 | prime | 159 |
| weierstrass | bn158 | prime | 158 |
| weierstrass | mnt1 | prime | 156 |
| weierstrass | secp128r1 | prime | 128 |
| weierstrass | secp128r2 | prime | 126 |
| weierstrass | secp112r1 | prime | 112 |
| weierstrass | secp112r2 | prime | 110 |

## Koblitz Form

| form | curve | field | n (bits) |
| --- | --- | --- | --- |
| koblitz | b571 | binary | 570 |
| koblitz | k571 | binary | 570 |
| koblitz | c2tnb431r1 | binary | 418 |
| koblitz | b409 | binary | 409 |
| koblitz | k409 | binary | 407 |
| koblitz | c2pnb368w1 | binary | 353 |
| koblitz | c2tnb359v1 | binary | 353 |
| koblitz | c2pnb304w1 | binary | 289 |
| koblitz | b283 | binary | 282 |
| koblitz | k283 | binary | 281 |
| koblitz | c2pnb272w1 | binary | 257 |
| koblitz | ansit239k1 | binary | 238 |
| koblitz | c2tnb239v1 | binary | 238 |
| koblitz | c2tnb239v2 | binary | 237 |
| koblitz | c2tnb239v3 | binary | 236 |
| koblitz | b233 | binary | 233 |
| koblitz | k233 | binary | 232 |
| koblitz | ansit193r1 | binary | 193 |
| koblitz | ansit193r2 | binary | 193 |
| koblitz | c2pnb208w1 | binary | 193 |
| koblitz | c2tnb191v1 | binary | 191 |
| koblitz | c2tnb191v2 | binary | 190 |
| koblitz | c2tnb191v3 | binary | 189 |
| koblitz | b163 | binary | 163 |
| koblitz | c2pnb163v1 | binary | 163 |
| koblitz | k163 | binary | 163 |
| koblitz | ansit163r1 | binary | 162 |
| koblitz | c2pnb163v2 | binary | 162 |
| koblitz | c2pnb163v3 | binary | 162 |
| koblitz | c2pnb176w1 | binary | 161 |
| koblitz | sect131r1 | binary | 131 |
| koblitz | sect131r2 | binary | 131 |
| koblitz | sect113r1 | binary | 113 |
| koblitz | sect113r2 | binary | 113 |
| koblitz | wap-wsg-idm-ecid-wtls1 | binary | 112 |

# Contributing

All PRs are more than welcome! If you are planning to contribute a large patch, please create an issue first to get any upfront questions or design decisions out of the way first.

You should be able run `make test` and `make lint` commands successfully before committing. Once a PR is created, GitHub test workflow will be run automatically and unit test results will be available in [GitHub actions](https://github.com/serengil/LightECC/actions/workflows/tests.yml) before approval.

# Support

There are many ways to support a project - starring⭐️ the GitHub repo is just one 🙏

You can also support this work on [Patreon](https://www.patreon.com/serengil?repo=lightecc), [GitHub Sponsors](https://github.com/sponsors/serengil) or [Buy Me a Coffee](https://buymeacoffee.com/serengil).

<a href="https://www.patreon.com/serengil?repo=lightecc">
<img src="https://raw.githubusercontent.com/serengil/LightPHE/master/icons/patreon.png" width="30%" height="30%">
</a>

<a href="https://buymeacoffee.com/serengil">
<img src="https://raw.githubusercontent.com/serengil/LightPHE/master/icons/bmc-button.png" width="25%" height="25%">
</a>

Also, your company's logo will be shown on README on GitHub if you become sponsor in gold, silver or bronze tiers.

# Citation

Please cite LightECC in your publications if it helps your research. Here is its BibTex entry:

```BibTeX
@misc{serengil2025lightecc
    author       = {Serengil, Sefik},
    title        = {LightECC: A Lightweight Elliptic Curve Cryptography Arithmetic Library for Python},
    year         = {2025},
    publisher    = {GitHub},
    howpublished = {\url{https://github.com/serengil/LightECC}},
}
```

# License

LightECC is licensed under the MIT License - see [`LICENSE`](https://github.com/serengil/LightECC/blob/master/LICENSE) for more details.

LightECC's [logo](https://thenounproject.com/icon/starfish-757257/) is designed by Identidea Portfolio.