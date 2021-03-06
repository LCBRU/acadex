# Acadex

Collate online information about academics

## Installation and Running

1. Download the code from github

```bash
git clone git@github.com:LCBRU/acadex.git
```

2. Install the requirements

Go to the `acadex` directory and type the command:

```bash
sudo apt install sqlite3
sudo apt-get install libldap2-dev
sudo apt-get install libsasl2-dev

pip install -r requirements.txt
```

3. Create the Environment File

Go to the `acadex` directory and copy the file `example.env` to `.env`.

Change the etails in the file to appropriate values.

4. Create the database using

Staying in the `acadex` directory and type the command:

```bash
./manage.py version_control
./manage.py upgrade
```

5. Run the application

Staying in the `acadex` directory and type the command:

```bash
./app.py
```

6. Start Celery Worker

This application uses Celery to run background tasks.
To start Celery run the following command from the `acadex`
directory:

```
celery -A celery_worker.celery worker
```

## Development

### Testing

To test the application, run the following command from the project folder:

```bash
pytest
```

### Database Schema Amendments

#### Create Migration

To create a migration run the command

```bash
./manage.py script "{Description of change}"
```

You will then need to change the newly-created script created in the
`migrations` directory to make the necessary upgrade and downgrade
changes.

#### Installation

To initialise the database run the commands:

```bash
manage.py version_control
manage.py upgrade
```

#### Upgrade

To upgrade the database to the current version, run the command:

```bash
manage.py upgrade
```
