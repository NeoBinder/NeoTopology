from mendeleev import element as mendeleev_element
from mendeleev import Element
import ttk


class UnknownElement:

    def __init__(self):
        self.symbol = "X"

    @property
    def mass(self):
        return 0.0


symbol_element_dict = {"X": UnknownElement()}


__all__ = ["Element"]


def element_from_symbol(symbol):
    symbol = ttk.data.NameUpperToSymbol.get(symbol, symbol)
    if symbol not in symbol_element_dict:
        try:
            element = mendeleev_element(symbol)
        except Exception:
            print(symbol)
            raise KeyError("symbol not found {}".format(symbol))
        symbol_element_dict[symbol] = element
    return symbol_element_dict[symbol]
