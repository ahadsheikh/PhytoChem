## PhytoChem Database

Source code for [phytochemdb.com](http://phytochemdb.com) built with Django and PostgreSQL database.

## Developer instructions

This project uses `conda` environment.
Make sure you have anaconda installed. Then you
have to create an environment with all necessary
packages. Use the following command to create
an environment from `environment.yml` file and activate
it-

```bash
conda env create -f environment.yml
conda activate PhytoChem
```

We have used PostgreSQL database. You need to setup a database first:

- Install postgresql in your machine: `sudo apt install postgresql` for Ubuntu and `sudo pacman -S postgresql` for Arch/Manjaro.
- Make sure the database service is running.
- Change to `postgres` user: `sudo su - postgres`
- Start `postgres` session: `psql`
- Now create a database: `CREATE DATABASE phytochem;`
- Create database user: `CREATE USER phytochemuser WITH PASSWORD '<Your Password>';`
- Modify some connection parameters:
  ```sql
  ALTER ROLE phytochemuser SET client_encoding TO 'utf8';
  ALTER ROLE phytochemuser SET default_transaction_isolation TO 'read committed';
  ALTER ROLE phytochemuser SET timezone TO 'UTC';
  ```
- Grant access rights to the database user: `GRANT ALL PRIVILEGES ON DATABASE phytochem TO phytochemuser;`
- Exit sql prompt: `\q`
- Exit user: `exit`

Create a `.env` file in the root directory specifying all the secrets as follows -

```
SECRET_KEY=1#e-7s-wm+)m!85^-b5-ymj+75wx*=df#5f6hi#_7(ji5n2nb(
DEBUG=True
DB_NAME=phytochem
DB_USER=phytochemuser
DB_PASS=<password>
DB_HOST=localhost
```
