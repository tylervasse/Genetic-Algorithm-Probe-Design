"""
IDT (Integrated DNA Technologies) API integration.

Handles authentication and melting-temperature lookups via the IDT OligoAnalyzer
REST API, plus the placeholder-strand design helpers that depend on those calls.
"""

import os

import requests
from Bio import Align, pairwise2

from dna_utils import rnatoDNA, revcomp

# Module-level token — refreshed automatically on 401/expiry.
access_token = None


def idtConnection():
    auth_url = "https://www.idtdna.com/Identityserver/connect/token"

    username = os.getenv("IDT_USERNAME")
    password = os.getenv("IDT_PASSWORD")
    client_id = os.getenv("IDT_CLIENT_ID")
    client_secret = os.getenv("IDT_CLIENT_SECRET")

    if not all([username, password, client_id, client_secret]):
        raise RuntimeError("Missing IDT API credentials in environment variables.")

    auth_data = {
        "grant_type": "password",
        "username": username,
        "password": password,
        "scope": "test",
    }

    auth_headers = {
        "Content-Type": "application/x-www-form-urlencoded",
    }

    auth_response = requests.post(
        auth_url,
        data=auth_data,
        headers=auth_headers,
        auth=(client_id, client_secret),
    )

    if auth_response.status_code == 200:
        return auth_response.json().get("access_token")

    raise RuntimeError(f"Authentication failed: {auth_response.status_code} {auth_response.text}")


def getCompMeltingTemp(sequence):
    """
    This function uses the IDT tools to calculate the melting temperature of a
    sequence with its complementary sequence. This function is used to find the
    maximum binding temperature of a miRNA to a placeholder strand and is also
    used to optimize the length of the complementary region of the placeholder
    strand to the probe.
    """
    global access_token

    try:
        analyze_url = "https://www.idtdna.com/restapi/v1/OligoAnalyzer/TmMisMatch"
        analyze_data = {
          "Settings": {
            "Sequence": sequence,
            "NaConc": 300,
            "MgConc": 0,
            "dNTPsConc": 0,
            "OligoConc": 0.2
          },
          "Sequence2": revcomp(sequence)
        }

        analyze_headers = {
                "Content-Type": "application/json",
                "Authorization": f"Bearer {access_token}",
            }

        analyze_response = requests.post(analyze_url, json=analyze_data, headers=analyze_headers)
        return analyze_response.json()["MinMeltTemp"]

    except:
        access_token = idtConnection()
        return getCompMeltingTemp(sequence)


def buildbasicplaceholder(miRNA):
    """
    The function builds the basic placeholder by providing the reverse complement
    of the miRNA.
    """
    miDNA = rnatoDNA(miRNA)
    placeholder = revcomp(miDNA)
    return placeholder


def getPlaceholderProbeMelt(MIRNA_PLACEHOLDER_MELT, probe_melting_temp):
    """
    If the automatic placeholder probe melting temperature is selected, this
    function will calculate an ideal placeholder probe melting temperature.
    It does this by averaging the melting temperature of the miRNA to the
    placeholder and the melting temperature of the probe. If this average is
    higher than 57C, the function returns 57C. If the probe melting temperature
    and the miRNA placeholder melting temperature are too far apart, the
    function simply returns the value of the probe melting temperature plus 5.
    """
    average = (MIRNA_PLACEHOLDER_MELT + probe_melting_temp) / 2

    if average >= 57.0 and probe_melting_temp < 50.0:
        return 57.0
    elif probe_melting_temp + 5 < MIRNA_PLACEHOLDER_MELT - 5:
        return probe_melting_temp + 5
    else:
        return average


def findPlaceholderSequences(miRNA, probe_sequence, probe_melting_temp, automatic_ph_mt, desired_ph_mt):
    """
    This function builds the placeholder sequences by calculating the melting
    temperature of the initial part of the reverse complement of the probe
    sequence. The function adds bases to the 3' end of the placeholder until
    the melting temperature is above the desired placeholder probe melting
    temperature. It then aligns this overlap portion of the placeholder with a
    basic design of a placeholder and adds missing bases to the 5' end until
    there is a 9 base toehold.
    """
    miRNA_placeholder_melt = getCompMeltingTemp(rnatoDNA(miRNA))
    rev_comp_probe_seq = revcomp(probe_sequence)
    base_index = 1
    overlap = rev_comp_probe_seq[:base_index]

    if automatic_ph_mt:
        placeholder_probe_melt = getPlaceholderProbeMelt(miRNA_placeholder_melt, probe_melting_temp)
    else:
        placeholder_probe_melt = desired_ph_mt

    while getCompMeltingTemp(overlap) < placeholder_probe_melt:
        base_index += 1
        overlap = rev_comp_probe_seq[:base_index]

    placeholder_probe_melting_temp = getCompMeltingTemp(overlap)

    matrix = Align.substitution_matrices.load("BLOSUM62")
    basic_placeholder = buildbasicplaceholder(miRNA)
    empty_probe = pairwise2.align.localds(overlap, basic_placeholder, matrix, -10, -0.5)[0][0]

    dash_end, dash_count = False, 0
    while not dash_end:
        if empty_probe[dash_count] != "-":
            dash_end = True
        else:
            dash_count += 1

    if dash_count == 10:
        placeholder = basic_placeholder[1:10] + overlap
        miRNA_placeholder_melt = getCompMeltingTemp(rnatoDNA(miRNA)[:-1])
    else:
        placeholder = basic_placeholder[:9] + overlap

    return placeholder, placeholder_probe_melting_temp, miRNA_placeholder_melt
