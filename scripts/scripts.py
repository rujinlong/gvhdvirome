#!/usr/bin/env python

def get_page_title(client, page_id, column="exp_sample_id"):
    """
    from notion_client import Client
    client = Client(auth=notion_token)
    """
    rst = client.pages.retrieve(page_id)
    return(rst['properties'][column]['title'][0]['plain_text'])