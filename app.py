import streamlit as st
from utils.cached_resources import load_collections
from utils.utilities import create_prompt


st.set_page_config(
    page_title="Auto Prompt",
    initial_sidebar_state="collapsed",
    page_icon="🤖",
    menu_items={"Report a Bug": "mailto:suh46@pitt.edu"},
)
vector_stores = load_collections()
search_query = st.text_area(
    "Question",
    placeholder="What's the best way to generate ALS treatment plan",
)
context_text = st.text_area(
    "Context",
    value="You are a ALZ clinician who knows how to generate ALS treatment plan",
)
col1, col2, col3 = st.columns(3)
with col1:
    prompt_size = st.number_input("prompt size", value=4096)
with col2:
    search_result_size = st.number_input("search result size", value=10)
with col3:
    collection_option = st.selectbox("From collection:", vector_stores.keys())

if len(search_query) > 0:
    search_results = vector_stores[collection_option].similarity_search_with_score(
        search_query, k=search_result_size
    )
    prompt = create_prompt(
        search_query,
        context_text,
        [(hit[0].metadata["source"], hit[0].page_content) for hit in search_results],
        prompt_size,
    )
    st.code(prompt, language=None)
