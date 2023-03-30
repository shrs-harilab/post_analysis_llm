import typing


def create_prompt(question_text:str, context_text:str, sources:typing.List[typing.Tuple[str, str]], prompt_length_limit:int = 1000):
    prompt_text = ""
    header_context = context_text + "\n\n"
    question_text = "\n" + question_text
    prompt_text += header_context

    words_left = prompt_length_limit - len(prompt_text) - len(question_text)

    for (source_id,text_segment) in sources:
        prompt_segment = """SOURCE:{0}\nTEXT:{1}\n\n""".format(source_id,text_segment)
        if len(prompt_segment) > words_left:
            break
        else:
            prompt_text += prompt_segment
    prompt_text += question_text

    return prompt_text


if __name__ == '__main__':
    create_prompt("I like this quesiton", [("source 1", "Acknowledgements We thank Dr. Rahim Moineddin for his assistance with statistical analysis.")])



