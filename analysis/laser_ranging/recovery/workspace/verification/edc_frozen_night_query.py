"""Retrieve the public EDC catalogue query for the frozen APOLLO night."""

from pathlib import Path
import requests
from bs4 import BeautifulSoup

root = Path("/home/dominic/Tudat/laser-ranging-work/research")
url = "https://edc.dgfi.tum.de/en/data/browser/"
session = requests.Session()
response = session.get(url, timeout=30)
response.raise_for_status()
soup = BeautifulSoup(response.text, "html.parser")
fields = {item["name"]: item.get("value", "") for item in soup.select("input[name]")}
for item in soup.select("select[name]"):
    option = item.select_one("option[selected]") or item.select_one("option")
    fields[item["name"]] = option.get("value", "") if option else ""
fields.update(
    query_data_type="NPT_v2",
    query_station="7045",
    query_start_data_date_operator="=",
    query_start_data_date_year="2017",
    query_start_data_date_month="3",
    query_start_data_date_day="9",
    query_show_out_of_date="yes",
    query_output_format="text",
)
response = session.post(url, data=fields, headers={"Referer": url}, timeout=20)
response.raise_for_status()
path = root / "EDC-APOL-20170309-v2-query.html"
path.write_bytes(response.content)
print(response.status_code, len(response.content), response.url)
print(BeautifulSoup(response.text, "html.parser").get_text(" ", strip=True)[-14000:])
