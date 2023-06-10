import typing
import tiktoken


def count_token(text: str):
    encoding = tiktoken.get_encoding("cl100k_base")
    encoded = encoding.encode(text)
    return len(encoded)


def create_prompt(
    question_text: str,
    context_text: str,
    sources: typing.List[typing.Tuple[str, str]],
    prompt_length_limit: int,
):
    prompt_text = ""
    header_context = context_text + "\n\n"
    question_text = "\n" + question_text
    prompt_text += header_context
    chars_left = (
        prompt_length_limit - count_token(prompt_text) - count_token(question_text)
    )
    for source_id, text_segment in sources:
        prompt_segment = """SOURCE:{0}\nTEXT:{1}\n\n""".format(source_id, text_segment)
        temp_tk = count_token(prompt_segment)
        if temp_tk > chars_left:
            continue
        else:
            prompt_text += prompt_segment
            chars_left -= temp_tk
    prompt_text += question_text

    return prompt_text
