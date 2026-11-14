"""Country-name normalisation to ISO 3166-1 alpha-3 codes (offline, no downloads).

Public repositories write countries many ways ("USA", "United States of America", "USA: Texas",
"Viet Nam"). Choropleth maps and cross-source joins need one code. Unknown names return ``None`` and
are *reported* by callers, never guessed.
"""

from __future__ import annotations

import re
import unicodedata

# ISO 3166-1 alpha-3 -> common English short name.
_ISO3 = """
ABW Aruba|AFG Afghanistan|AGO Angola|AIA Anguilla|ALA Åland Islands|ALB Albania|AND Andorra|ARE United Arab Emirates|
ARG Argentina|ARM Armenia|ASM American Samoa|ATA Antarctica|ATF French Southern Territories|ATG Antigua and Barbuda|
AUS Australia|AUT Austria|AZE Azerbaijan|BDI Burundi|BEL Belgium|BEN Benin|BES Bonaire, Sint Eustatius and Saba|
BFA Burkina Faso|BGD Bangladesh|BGR Bulgaria|BHR Bahrain|BHS Bahamas|BIH Bosnia and Herzegovina|BLM Saint Barthélemy|
BLR Belarus|BLZ Belize|BMU Bermuda|BOL Bolivia|BRA Brazil|BRB Barbados|BRN Brunei|BTN Bhutan|BVT Bouvet Island|
BWA Botswana|CAF Central African Republic|CAN Canada|CCK Cocos (Keeling) Islands|CHE Switzerland|CHL Chile|CHN China|
CIV Côte d'Ivoire|CMR Cameroon|COD Democratic Republic of the Congo|COG Republic of the Congo|COK Cook Islands|
COL Colombia|COM Comoros|CPV Cabo Verde|CRI Costa Rica|CUB Cuba|CUW Curaçao|CXR Christmas Island|CYM Cayman Islands|
CYP Cyprus|CZE Czechia|DEU Germany|DJI Djibouti|DMA Dominica|DNK Denmark|DOM Dominican Republic|DZA Algeria|
ECU Ecuador|EGY Egypt|ERI Eritrea|ESH Western Sahara|ESP Spain|EST Estonia|ETH Ethiopia|FIN Finland|FJI Fiji|
FLK Falkland Islands|FRA France|FRO Faroe Islands|FSM Micronesia|GAB Gabon|GBR United Kingdom|GEO Georgia|
GGY Guernsey|GHA Ghana|GIB Gibraltar|GIN Guinea|GLP Guadeloupe|GMB Gambia|GNB Guinea-Bissau|GNQ Equatorial Guinea|
GRC Greece|GRD Grenada|GRL Greenland|GTM Guatemala|GUF French Guiana|GUM Guam|GUY Guyana|HKG Hong Kong|
HMD Heard Island and McDonald Islands|HND Honduras|HRV Croatia|HTI Haiti|HUN Hungary|IDN Indonesia|IMN Isle of Man|
IND India|IOT British Indian Ocean Territory|IRL Ireland|IRN Iran|IRQ Iraq|ISL Iceland|ISR Israel|ITA Italy|
JAM Jamaica|JEY Jersey|JOR Jordan|JPN Japan|KAZ Kazakhstan|KEN Kenya|KGZ Kyrgyzstan|KHM Cambodia|KIR Kiribati|
KNA Saint Kitts and Nevis|KOR South Korea|KWT Kuwait|LAO Laos|LBN Lebanon|LBR Liberia|LBY Libya|LCA Saint Lucia|
LIE Liechtenstein|LKA Sri Lanka|LSO Lesotho|LTU Lithuania|LUX Luxembourg|LVA Latvia|MAC Macao|MAF Saint Martin|
MAR Morocco|MCO Monaco|MDA Moldova|MDG Madagascar|MDV Maldives|MEX Mexico|MHL Marshall Islands|MKD North Macedonia|
MLI Mali|MLT Malta|MMR Myanmar|MNE Montenegro|MNG Mongolia|MNP Northern Mariana Islands|MOZ Mozambique|
MRT Mauritania|MSR Montserrat|MTQ Martinique|MUS Mauritius|MWI Malawi|MYS Malaysia|MYT Mayotte|NAM Namibia|
NCL New Caledonia|NER Niger|NFK Norfolk Island|NGA Nigeria|NIC Nicaragua|NIU Niue|NLD Netherlands|NOR Norway|
NPL Nepal|NRU Nauru|NZL New Zealand|OMN Oman|PAK Pakistan|PAN Panama|PCN Pitcairn|PER Peru|PHL Philippines|
PLW Palau|PNG Papua New Guinea|POL Poland|PRI Puerto Rico|PRK North Korea|PRT Portugal|PRY Paraguay|PSE Palestine|
PYF French Polynesia|QAT Qatar|REU Réunion|ROU Romania|RUS Russia|RWA Rwanda|SAU Saudi Arabia|SDN Sudan|
SEN Senegal|SGP Singapore|SGS South Georgia and the South Sandwich Islands|SHN Saint Helena|SJM Svalbard and Jan Mayen|
SLB Solomon Islands|SLE Sierra Leone|SLV El Salvador|SMR San Marino|SOM Somalia|SPM Saint Pierre and Miquelon|
SRB Serbia|SSD South Sudan|STP Sao Tome and Principe|SUR Suriname|SVK Slovakia|SVN Slovenia|SWE Sweden|SWZ Eswatini|
SXM Sint Maarten|SYC Seychelles|SYR Syria|TCA Turks and Caicos Islands|TCD Chad|TGO Togo|THA Thailand|
TJK Tajikistan|TKL Tokelau|TKM Turkmenistan|TLS Timor-Leste|TON Tonga|TTO Trinidad and Tobago|TUN Tunisia|
TUR Türkiye|TUV Tuvalu|TWN Taiwan|TZA Tanzania|UGA Uganda|UKR Ukraine|UMI United States Minor Outlying Islands|
URY Uruguay|USA United States|UZB Uzbekistan|VAT Holy See|VCT Saint Vincent and the Grenadines|VEN Venezuela|
VGB British Virgin Islands|VIR U.S. Virgin Islands|VNM Vietnam|VUT Vanuatu|WLF Wallis and Futuna|WSM Samoa|
XKX Kosovo|YEM Yemen|ZAF South Africa|ZMB Zambia|ZWE Zimbabwe
"""

# Extra spellings seen in public repositories (key after normalisation -> ISO3).
_ALIASES = {
    "usa": "USA",
    "us": "USA",
    "u s a": "USA",
    "united states of america": "USA",
    "america": "USA",
    "uk": "GBR",
    "great britain": "GBR",
    "england": "GBR",
    "scotland": "GBR",
    "wales": "GBR",
    "northern ireland": "GBR",
    "united kingdom of great britain and northern ireland": "GBR",
    "russian federation": "RUS",
    "viet nam": "VNM",
    "turkey": "TUR",
    "turkiye": "TUR",
    "korea": "KOR",
    "republic of korea": "KOR",
    "korea south": "KOR",
    "south korea": "KOR",
    "korea north": "PRK",
    "democratic people s republic of korea": "PRK",
    "czech republic": "CZE",
    "ivory coast": "CIV",
    "cote d ivoire": "CIV",
    "cape verde": "CPV",
    "swaziland": "SWZ",
    "macedonia": "MKD",
    "the netherlands": "NLD",
    "holland": "NLD",
    "burma": "MMR",
    "east timor": "TLS",
    "laos pdr": "LAO",
    "lao people s democratic republic": "LAO",
    "lao pdr": "LAO",
    "syrian arab republic": "SYR",
    "iran islamic republic of": "IRN",
    "islamic republic of iran": "IRN",
    "bolivia plurinational state of": "BOL",
    "venezuela bolivarian republic of": "VEN",
    "tanzania united republic of": "TZA",
    "united republic of tanzania": "TZA",
    "moldova republic of": "MDA",
    "republic of moldova": "MDA",
    "palestinian territory": "PSE",
    "state of palestine": "PSE",
    "west bank": "PSE",
    "gaza strip": "PSE",
    "democratic republic of congo": "COD",
    "congo democratic republic": "COD",
    "drc": "COD",
    "dr congo": "COD",
    "congo kinshasa": "COD",
    "congo brazzaville": "COG",
    "congo": "COG",
    "the gambia": "GMB",
    "the bahamas": "BHS",
    "brunei darussalam": "BRN",
    "micronesia federated states of": "FSM",
    "hong kong sar": "HKG",
    "macau": "MAC",
    "reunion": "REU",
    "curacao": "CUW",
    "saint helena ascension and tristan da cunha": "SHN",
    "sao tome and principe": "STP",
    "timor leste": "TLS",
    "vatican": "VAT",
    "vatican city": "VAT",
    "falkland islands malvinas": "FLK",
    "aland islands": "ALA",
    "korea republic of": "KOR",
    "korea democratic people s republic of": "PRK",
    "congo the democratic republic of the": "COD",
    "congo democratic republic of the": "COD",
    "taiwan province of china": "TWN",
    "virgin islands british": "VGB",
    "virgin islands u s": "VIR",
    "holy see vatican city state": "VAT",
    "sudan south": "SSD",
    "slovak republic": "SVK",
    "kyrgyz republic": "KGZ",
    "republic of congo": "COG",
    "svalbard": "SJM",
    "gambia the": "GMB",
    "hong kong sar china": "HKG",
    "libyan arab jamahiriya": "LBY",
    "macao sar china": "MAC",
    "russian federation the": "RUS",
    "turkey republic of": "TUR",
    "eswatini kingdom of": "SWZ",
}


def _key(text: object) -> str:
    """Accent-, case- and punctuation-insensitive lookup key."""
    value = unicodedata.normalize("NFKD", str(text))
    value = "".join(ch for ch in value if not unicodedata.combining(ch))
    value = value.lower().replace("&", " and ")
    value = re.sub(r"[^a-z0-9]+", " ", value).strip()
    value = re.sub(r"\bst\b", "saint", value)
    return re.sub(r"^the ", "", value)


def _build() -> tuple[dict[str, str], dict[str, str]]:
    names: dict[str, str] = {}
    lookup: dict[str, str] = {}
    for entry in re.split(r"\|\s*", _ISO3.replace("\n", " ").strip()):
        entry = entry.strip()
        if not entry:
            continue
        code, name = entry.split(" ", 1)
        names[code] = name.strip()
        lookup[_key(name)] = code
        lookup[code.lower()] = code
    for alias, code in _ALIASES.items():
        lookup[_key(alias)] = code
    return names, lookup


ISO3_NAMES, _LOOKUP = _build()


def country_to_iso3(value: object) -> str | None:
    """ISO 3166-1 alpha-3 for a country string, or ``None`` when it cannot be resolved.

    Handles the repository convention ``"Country: region, locality"`` by using the part before the
    first colon, and ignores a trailing parenthetical.
    """
    if value is None:
        return None
    try:
        if value != value:  # NaN
            return None
    except (TypeError, ValueError):
        return None
    text = str(value).strip()
    if not text:
        return None
    head = text.split(":", 1)[0]
    candidates = [head, re.sub(r"\(.*?\)", "", head)]
    for candidate in candidates:
        key = _key(candidate)
        if key in _LOOKUP:
            return _LOOKUP[key]
    return None


def iso3_to_name(code: str) -> str | None:
    return ISO3_NAMES.get(str(code).upper())
