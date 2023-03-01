# ALS prompt enigneer



## Getting start
Create an virtual environment for python 3.8.8

Install dependencies
```shell
pip install -r requirements.txt
```

create a env file and update the environment variables. (Just change it to your Mongodb's username and password)
```shell
cp .env.example .env
```

create milvus database from docker-compose
```shell
docker compose up
```

## Prepare data
* run ```load_data.py``` to download the data from mongodb to a local pickle file
* run ```load_to_milvus.py``` to load the downloaded data into milvus
* run ```python -m streamlit run analyze_als.py``` and it will open the semantic search interface