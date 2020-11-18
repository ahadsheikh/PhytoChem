## PhytoChem Dataset
Store and show chemical dataset using Django

## Installation instruction
This project uses `conda` environment.
Make sure you have anaconda installed. Then you
have to create an environment with all necessary
libraries. Use the following command to create
an environment from `environment.yml` file and activate
it-

```
conda env create -f environment.yml
conda activate PhytoChem
```

Can be used with `heroku` with this buildpack: [https://github.com/pl31/heroku-buildpack-conda](https://github.com/pl31/heroku-buildpack-conda), but our project is too large after compression, so bad luck :(
