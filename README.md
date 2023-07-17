# Fragment hit follow chemistry collection

> :construction: This repository is under construction :construction:

A collection of notebooks and scripts for the prediction of follow-up compounds 
I used for follow-up chemistry for XChem targets,
_ie._ small fragments identified crystallographically.

It is not intended as an exhaustive repository of all possible methods.

For now see: https://github.com/matteoferla/EV-D68-3C-protease

## XChem libraries

[combined XChem libraries](combined-XChem-libraries.csv) is a CSV file of the compounds in the XChem libraries.

## Fragment theft

...

## Fragment extraction

...
## Fragment PLIP

...
## Fragmenstein

> For more [:octocat::Fragmenstein](https://github.com/matteoferla/Fragmenstein)

...
## Fragment network

...
## Arthorian quest

> See [:octocat::Arthorian quest](https://github.com/matteoferla/arthorian-quest)
> (:construction: under construction :construction:)

...

## Filtering

...


## CoPriNet

> See [:octocat::CoPriNet](https://github.com/oxpig/CoPriNet)

Predict the price for a compound.

...

## Fragalysis

Typical folder format.

### Fauxalysis

### Fragalysis Upload

This is just to prep.
The upload form is at https://fragalysis.diamond.ac.uk/viewer/upload_cset/
> :warning: The upload code exists in Fragalysis-API old version. Hunt it down.

```python
from followup.prep_fragalysis import prep, generate_header
header = generate_header(method='Foo',
                         ref_url='https://www.example.com',
                         submitter_name='unknown',
                         submitter_email='a@b.c',
                         submitter_institution='Nowehere',
                         generation_date='2012-12-12',
                         smiles='CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                         extras={'∆∆G': 'W. Josiah Gibbs'})
prep(df, header, 'mol', 'id', 'foo.sdf', 'x1234', 'x1234', ['∆∆G'])
```

## Michelanglo

Upload to Fragalysis is not always okay.

## Minor

### Zinc data

`ZincInformer` in [zinc_data.py](followup/zinc_data.py) is a class to get data from a series of ZINC IDs.

