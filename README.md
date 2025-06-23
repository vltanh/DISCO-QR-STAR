# DISCO-QR-STAR
UIUC CS 581 Final Project

## Usage

### Run

### Visualization

```
python check_missing.py --input output/trees/ --reference data/trees/
```

To aggregate the results, run
```
python agg_result.py --input output/trees/ --output ncd.csv --reference data/trees/
```

To visualize, run
```
python plot_result.py --input ncd.csv --output plots/ --species-tree trues
```