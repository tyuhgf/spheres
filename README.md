#
Tool for working with sphere triangulations.
Can be used to calculate first Pontryagin class of a combinatorial manifold.

## Testing
For all tests:
```
python -m pytest tests --cov=spheres
```

For tests not involving GAP execution:
```
python -m pytest tests --cov=spheres -m "not gap"
```

## Preparation

1. Install GAP (https://www.gap-system.org/Download/)
2. ```
    conda install -c conda-forge sage  
    pip install -r requirements.txt
    ```
3. copy `settings/local.py.default` to `settings/local.py`, set the valid path to GAP executable

