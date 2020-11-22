## PhytoChem Database
Source code for [phytochemdb.com](http://phytochemdb.com) built with Django

## Developer instructions
This project uses `conda` environment.
Make sure you have anaconda installed. Then you
have to create an environment with all necessary
packages. Use the following command to create
an environment from `environment.yml` file and activate
it-

```
conda env create -f environment.yml
conda activate PhytoChem
```

You also need to create a `.env` file in the parent directory like this, the secret key is different on the host server -

```
SECRET_KEY=1#e-7s-wm+)m!85^-b5-ymj+75wx*=df#5f6hi#_7(ji5n2nb(
DEBUG=True
```
