# prompt enigneer for ALS online post analysis



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
  
## Start semantic search interface
* run ```python -m streamlit run analyze_als.py``` and it will open the semantic search interface


# knowledge base
## Prepare data
* run ```pip install -r requirements.txt``` to install new dependecies
* create a folder called ```pdfs``` inside knowledge_base, put the pdf files you want to index inside it
* run ```knowledge_base/load_to_milvus.py``` to extract pdf files from ```pdfs``` and load to milvus database

## Start auto prompt creation UI
* run ```python -m streamlit run knowledge_base/auto_generate_prompt.py```
```
