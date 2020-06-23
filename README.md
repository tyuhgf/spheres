#
Tool for working with sphere triangulations based on Sage, Bistellar, ...

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
```
conda install -c conda-forge sage  
pip install -r requirements.txt
```