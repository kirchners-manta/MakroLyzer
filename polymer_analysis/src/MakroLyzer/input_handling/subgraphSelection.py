# for handling regular expressions
import re

def parse_selection_list(selection):
    """
    Parse a selection list into normalized selection dicts.

    Selection tokens:
      1) Exact composition: contains at least one digit (e.g. C2H6, CH4, C1H4).
      2) Element-only: no digits (e.g. CHO) -> match any counts with exactly these elements.

    Args:
        selection (str | list[str]): Selection tokens or a comma-separated string.

    Returns:
        list[dict]: Normalized selections with keys: mode, elements, counts (for exact).
    """
    if selection is None:
        return []
    
    # check if we got only one or multiple tokens
    if isinstance(selection, str):
        tokens = [selection]
    else:
        tokens = list(selection)
        
    normalized = []
    for token in tokens:
        if token is None:
            continue
        norm_token = _parse_selection_token(token)
        normalized.append(norm_token)
    return normalized

def _parse_selection_token(token):
    token = token.strip()
    if not token:
        raise ValueError("Empty selection token.")
    
    # check if we got option 1 or option 2
    has_digit = any(pos.isdigit() for pos in token)
    
    elements = []
    counts = {} # dictionary
    
    i=0
    while i < len(token):
        # check if the element is availible in MakroLyzer.dictionaries
        # re checks if a pattern starting at position i in the token matches on of the elements
        # long element names are in the front of the check to prevent for example Cl to be identified
        # as C
        match = re.match(r'(Li|Na|K|Mg|Ca|Zn|Cl|Br|Ag|Au|N|O|B|F|P|S|C|H|I|D|X)', token[i:])
        if not match:
            raise ValueError(f"Invalid selection token: {token}")
    
        element = match.group(0)
        i+=len(element) # jump forward in the token
        
        # read the single digits between the letters and combine them to obtain the count
        number = []
        while i < len(token) and token[i].isdigit():
            number.append(token[i])
            i+=1 # jump forward in the token
        count = int(''.join(number)) if number else 1
        
        elements.append(element)
        counts[element] = counts.get(element, 0) + count # we add the count in case the element appears multiple times
        
    if has_digit:
        return {
            "mode": "exact",
            "elements": set(counts.keys()),
            "counts": counts,
            "raw": token,
        }

    return {
        "mode": "elements",
        "elements": set(elements),
        "counts": None,
        "raw": token,
    } 