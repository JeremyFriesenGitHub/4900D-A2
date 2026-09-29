# Quantitative results

2000 reconstruction samples and 20000 independent test points per mesh; meshes normalized to [-1, 1]; marching cubes on a 64^3 grid over [-1.5, 1.5]^3.
Point-to-plane: phi(r) = 1/(r^6 + eps^6) with eps = 0.5 h (h = mean kNN spacing of the samples). Winding number: k = 6.

## 1. Error |Phi(x)| at the test points, clean samples (required comparison)

Units differ: point-to-plane Phi is a length, 0.5 - w is a fraction of the full solid angle.

| mesh | method | max \|Phi\| | mean \|Phi\| | median \|Phi\| |
|---|---|---:|---:|---:|
| sphere | winding number | 4.6560 | 0.0272 | 0.0156 |
| sphere | point-to-plane | 0.0313 | 0.0017 | 0.0007 |
| cube | winding number | 10.3887 | 0.0393 | 0.0123 |
| cube | point-to-plane | 0.0998 | 0.0025 | 0.0001 |
| bunny | winding number | 125.9214 | 0.0784 | 0.0378 |
| bunny | point-to-plane | 0.0566 | 0.0032 | 0.0015 |
| camel | winding number | 139.2207 | 0.0892 | 0.0384 |
| camel | point-to-plane | 0.0555 | 0.0037 | 0.0018 |
| kid | winding number | 43.8899 | 0.0818 | 0.0427 |
| kid | point-to-plane | 0.0291 | 0.0023 | 0.0011 |

## 2. Distances measured on the reconstructed meshes, clean samples

test->recon: test point to reconstructed surface. recon->truth: reconstructed vertex to true surface (large values = spurious surfaces).

| mesh | method | test->recon max | test->recon mean | recon->truth max | recon->truth mean |
|---|---|---:|---:|---:|---:|
| sphere | winding number | 0.0308 | 0.0040 | 0.0410 | 0.0045 |
| sphere | point-to-plane | 0.0294 | 0.0015 | 0.0315 | 0.0017 |
| cube | winding number | 0.0855 | 0.0047 | 0.0592 | 0.0055 |
| cube | point-to-plane | 0.0871 | 0.0024 | 0.1709 | 0.0077 |
| bunny | winding number | 0.0611 | 0.0055 | 0.0980 | 0.0063 |
| bunny | point-to-plane | 0.0450 | 0.0026 | 0.2060 | 0.0070 |
| camel | winding number | 0.0882 | 0.0077 | 0.0403 | 0.0060 |
| camel | point-to-plane | 0.0564 | 0.0033 | 0.3455 | 0.0230 |
| kid | winding number | 0.0847 | 0.0074 | 0.0348 | 0.0063 |
| kid | point-to-plane | 0.0500 | 0.0027 | 1.2420 | 0.2782 |

## 3. Reference variants, clean samples

| mesh | method | max \|Phi\| | mean \|Phi\| | median \|Phi\| | test->recon max | test->recon mean | recon->truth max | recon->truth mean |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| sphere | point-to-plane, handout weight | 0.4207 | 0.2063 | 0.2115 | no surface | no surface | no surface | no surface |
| sphere | point-to-plane, handout weight, 16 nearest | 0.0258 | 0.0040 | 0.0032 | 0.0251 | 0.0038 | 0.0263 | 0.0041 |
| sphere | point-to-plane, squared weight | 0.0763 | 0.0049 | 0.0027 | 0.0827 | 0.0051 | 0.1014 | 0.0057 |
| sphere | point-to-plane, exponential weight | 0.0254 | 0.0033 | 0.0029 | 0.0198 | 0.0029 | 0.0218 | 0.0032 |
| sphere | winding number, regularized | 0.0823 | 0.0170 | 0.0138 | 0.0161 | 0.0021 | 0.0130 | 0.0021 |
| cube | point-to-plane, handout weight | 0.4256 | 0.2309 | 0.2321 | no surface | no surface | no surface | no surface |
| cube | point-to-plane, handout weight, 16 nearest | 0.0926 | 0.0059 | 0.0000 | 0.1149 | 0.0070 | 0.1426 | 0.0112 |
| cube | point-to-plane, squared weight | 0.0976 | 0.0072 | 0.0033 | 0.1030 | 0.0082 | 0.1743 | 0.0140 |
| cube | point-to-plane, exponential weight | 0.0965 | 0.0053 | 0.0001 | 0.0870 | 0.0055 | 0.2854 | 0.0101 |
| cube | winding number, regularized | 0.3741 | 0.0365 | 0.0145 | 0.0943 | 0.0050 | 0.0625 | 0.0041 |
| bunny | point-to-plane, handout weight | 0.2641 | 0.1242 | 0.1262 | no surface | no surface | no surface | no surface |
| bunny | point-to-plane, handout weight, 16 nearest | 0.0546 | 0.0062 | 0.0042 | 0.0779 | 0.0067 | 0.1890 | 0.0125 |
| bunny | point-to-plane, squared weight | 0.0539 | 0.0053 | 0.0029 | 0.0557 | 0.0058 | 0.3414 | 0.0168 |
| bunny | point-to-plane, exponential weight | 0.0521 | 0.0060 | 0.0043 | 0.0487 | 0.0052 | 1.2036 | 0.1192 |
| bunny | winding number, regularized | 0.4016 | 0.0542 | 0.0386 | 0.1021 | 0.0062 | 0.0972 | 0.0053 |
| camel | point-to-plane, handout weight | 0.1550 | 0.0627 | 0.0594 | 3.5896 | 2.4858 | 1.5944 | 1.5724 |
| camel | point-to-plane, handout weight, 16 nearest | 0.0669 | 0.0072 | 0.0049 | 0.0912 | 0.0086 | 1.3539 | 0.1505 |
| camel | point-to-plane, squared weight | 0.0496 | 0.0057 | 0.0034 | 0.0646 | 0.0069 | 1.8434 | 0.3582 |
| camel | point-to-plane, exponential weight | 0.0594 | 0.0069 | 0.0048 | 0.0624 | 0.0057 | 1.5447 | 0.3472 |
| camel | winding number, regularized | 0.4566 | 0.0815 | 0.0521 | 0.1265 | 0.0115 | 0.0369 | 0.0050 |
| kid | point-to-plane, handout weight | 0.0845 | 0.0401 | 0.0397 | 2.4694 | 1.6242 | 1.6615 | 1.1557 |
| kid | point-to-plane, handout weight, 16 nearest | 0.0375 | 0.0049 | 0.0035 | 0.0501 | 0.0043 | 1.9736 | 0.5770 |
| kid | point-to-plane, squared weight | 0.0348 | 0.0038 | 0.0022 | 0.0496 | 0.0053 | 2.0015 | 0.6199 |
| kid | point-to-plane, exponential weight | 0.0347 | 0.0043 | 0.0034 | 0.0301 | 0.0026 | 1.8636 | 0.6721 |
| kid | winding number, regularized | 0.5448 | 0.0702 | 0.0573 | 0.1112 | 0.0084 | 0.0239 | 0.0044 |

## 4. Noisy samples (sigma = 0.02); test points stay on the clean surface

| mesh | method | max \|Phi\| | mean \|Phi\| | median \|Phi\| | test->recon max | test->recon mean | recon->truth max | recon->truth mean |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| sphere | winding number | 383.1845 | 0.3588 | 0.1385 | 0.0527 | 0.0082 | 0.0773 | 0.0105 |
| sphere | point-to-plane | 0.0497 | 0.0098 | 0.0081 | 0.0517 | 0.0094 | 0.0647 | 0.0120 |
| cube | winding number | 370.7716 | 0.3923 | 0.1182 | 0.0752 | 0.0087 | 0.1019 | 0.0119 |
| cube | point-to-plane | 0.0721 | 0.0116 | 0.0094 | 0.0751 | 0.0113 | 0.1705 | 0.0181 |
| bunny | winding number | 821.5820 | 0.4112 | 0.1708 | 0.0737 | 0.0090 | 0.1007 | 0.0122 |
| bunny | point-to-plane | 0.0553 | 0.0101 | 0.0085 | 0.0623 | 0.0098 | 0.1771 | 0.0170 |
| camel | winding number | 66.2386 | 0.3468 | 0.1956 | 0.1016 | 0.0111 | 0.0732 | 0.0126 |
| camel | point-to-plane | 0.0539 | 0.0091 | 0.0074 | 0.0601 | 0.0095 | 0.4409 | 0.0448 |
| kid | winding number | 116.4863 | 0.3764 | 0.2369 | 0.0774 | 0.0112 | 0.0955 | 0.0144 |
| kid | point-to-plane | 0.0349 | 0.0069 | 0.0057 | 0.0406 | 0.0083 | 1.1499 | 0.1695 |
