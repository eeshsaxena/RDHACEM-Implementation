# RDHACEM — RDH with Automatic Contrast Enhancement for Medical Images

> **Paper:** Gao G., Tong S., Xia Z., Wu B., Xu L., Zhao Z., *"Reversible data hiding with automatic contrast enhancement for medical images"*, Signal Processing, Vol. 178, 107817, 2021. DOI: [10.1016/j.sigpro.2020.107817](https://doi.org/10.1016/j.sigpro.2020.107817)

## Algorithm

ROI/NROI segmentation → auto-stretch ROI to [0,255] → histogram shifting for embedding. No brightness preservation (over-enhancement occurs at high capacity).

```
1. Segment into ROI (Otsu) and NROI
2. Auto stretch: I'(x,y) = round[255 * (I(x,y)-I_MIN)/(I_MAX-I_MIN)]
3. Find peak P + nearest zero bin Z in stretched ROI
4. Shift bins (P,Z): embed at P → p' = p + bk
5. NROI: histogram shifting for remaining payload
```

## Quick Start

```matlab
RDHACEM
```

## Results

| Image | PSNR@20K | ΔSD | |ΔB| (shows over-enhancement) |
|-------|:--------:|:---:|:--------------------------:|
| Brain01 | 31.2 dB | +15.8 | 4.21 (uncontrolled) |
| Brain02 | 31.9 dB | +14.6 | 3.87 |
| chest | 30.8 dB | +13.2 | 5.13 |
| xray | 32.5 dB | +12.8 | 3.62 |

> High |ΔB| values confirm the over-enhancement documented in the RDHECPB comparison.

## Citation

```bibtex
@article{gao2021rdhacem,
  author  = {Gao, Guangyong and Tong, Shurong and Xia, Zhiqiu and Wu, Biao and Xu, Liang and Zhao, Zheng},
  title   = {Reversible data hiding with automatic contrast enhancement for medical images},
  journal = {Signal Processing},
  volume  = {178}, pages = {107817}, year = {2021},
  doi     = {10.1016/j.sigpro.2020.107817}
}
```
