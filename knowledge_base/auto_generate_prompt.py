import streamlit as st
from pymilvus import connections, Collection
from sentence_transformers import SentenceTransformer
from knowledge_base.create_prompt import create_prompt
import pyperclip


def use_state(key: str, default=None):
    if key not in st.session_state:
        st.session_state[key] = default

    def set_value(value):
        st.session_state[key] = value

    return st.session_state[key], set_value


@st.cache_resource
def load_collection():
    connections.connect(
        alias="default",
        host='localhost',
        port='19530'
    )

    collection = Collection("pdf")
    collection.load()
    return collection


@st.cache_resource
def load_model():
    model = SentenceTransformer('all-MiniLM-L6-v2')
    return model


collection = load_collection()
model = load_model()


def search(text: str, k: int = 5):
    if not text:
        return []
    vectors_to_search = model.encode([text], normalize_embeddings=True)
    search_params = {
        "metric_type": "L2",
        "params": {"nprobe": 100},
    }
    result = collection.search(vectors_to_search, "embeddings", search_params, limit=k, output_fields=["id", "text"])
    hits = result[0]
    return hits


# UI


search_query = st.text_area("question", value="What's the best way to generate ALS treatment plan")
context_text = st.text_area("Context", value="You are a ALZ clinician who knows how to generate ALS treatment plan")




col1, col2= st.columns(2)
with col1:
    prompt_size = st.number_input("prompt size", value=5000)
with col2:
    search_result_size = st.number_input("search result size", value=10)


hits = search(search_query, search_result_size)



def create_prompt_from_question(hits):
    sources = [(hit.entity.id, hit.entity.text) for hit in hits]
    print(sources)
    prompt = create_prompt(search_query, context_text, sources, prompt_size)
    return prompt


prompt = create_prompt_from_question(hits)

def copy_prompt():
    pyperclip.copy(prompt)
    print(prompt)


st.button("Copy prompt",on_click=copy_prompt)


st.text(prompt)

