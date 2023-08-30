# HoMi
Host-microbiome dual transcriptome pipeline

## Installation
```
git clone git@github.com:sterrettJD/HoMi.git
cd HoMi
pip install -e .
```

## Usage
```
HoMi.py <config_file>
```

### Running with example dataset
```
# Create the mock community (too big for github)
python create_mock_community.py

# Use the example config and example metadata provided
HoMi.py tests/example_config.yaml
```