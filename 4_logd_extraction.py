import pandas as pd
import requests

df = pd.read_csv("unique_smiles.csv")

# API url for ADMETlab3, mimics browser request to get logD7.4 values for each SMILES
url = "https://admetlab3.scbdd.com/server/evaluationCal"
headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:143.0) Gecko/20100101 Firefox/143.0",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.5",
    "Content-Type": "application/x-www-form-urlencoded",
    "Sec-GPC": "1",
    "Upgrade-Insecure-Requests": "1",
    "Sec-Fetch-Dest": "document",
    "Sec-Fetch-Mode": "navigate",
    "Sec-Fetch-Site": "same-origin",
    "Sec-Fetch-User": "?1",
    "Referer": "https://admetlab3.scbdd.com/server/evaluation",
    "Priority": "u=0, i",
}
data = {
    "csrfmiddlewaretoken": "2dbRH9yczyItLmYWoPUAKeSIpbqhjmI375HPKTIwJvNR7bU0Fwc53Dq5khs0XtPk",
    "method": 1,
}

results = pd.DataFrame(columns=["smiles", "logD7.4"])

for index, row in df.iterrows():
    smiles = row["smiles"]
    print(smiles)
    data["smiles"] = smiles

    # Ignore expired SSL certificate
    response = requests.post(url, headers=headers, data=data, verify=False)

    logd_index = response.text.find("logD7.4")
    td_index = response.text.find("<td", logd_index)
    td_index = response.text.find(">", td_index)
    end_td_index = response.text.find("</td>", td_index)

    logd = response.text[td_index + 1 : end_td_index]
    results = pd.concat(
        [results, pd.DataFrame({"smiles": [smiles], "logD7.4": [logd]})],
        ignore_index=True,
    )
    print(logd)

results.to_csv("results_logd.csv", index=False)
print("Results saved")
