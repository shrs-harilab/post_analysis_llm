# Prompt Enigneer and Semantic Search for Multiple Knowledge or Forums

## Getting Start
Create an virtual environment for python

Install dependencies in the virtual environment
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

## Forums
### Prepare Data for ALZ
* Run ```forums/alz/load_data.py``` to download the data from mongodb to a local pickle file
* Run ```forums/alz/load_to_milvus.py``` to load the downloaded data into milvus
## Knowledge base
### Prepare data for PubMed (PMC)
* Run ```knowledge_base/NCBI/Crawler.py``` to download (scrape) pdfs from PMC to ```knowledge_base/NCBI/pdfs/```
* Run ```knowledge_base/NCBI/load_to_milvus.py``` to load the pdfs into milvus
## Start the Auto Prompt Creation UI
* run ```python -m streamlit run app.py```

