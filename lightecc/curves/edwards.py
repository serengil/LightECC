# project dependencies
from lightecc.interfaces.form import TwistedEdwardsInterface
from lightecc.commons.logger import Logger

logger = Logger(module="lightecc/curves/edwards.py")

DEFAULT_CURVE = "ed25519"


# pylint: disable=too-few-public-methods
class Ed25519(TwistedEdwardsInterface):
    p = pow(2, 255) - 19
    a = -1
    d = (-121665 * pow(121666, -1, p)) % p

    u = 9
    g_y = ((u - 1) * pow(u + 1, -1, p)) % p
    g_x = 15112221349535400772501151409588531511454012693041857206046113283949847762202
    G = (g_x, g_y)

    n = 0x1000000000000000000000000000000014DEF9DEA2F79CD65812631A5CF5D3ED


class Ed448(TwistedEdwardsInterface):
    p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    a = 1
    d = 0xD78B4BDC7F0DAF19F24F38C29373A2CCAD46157242A50F37809B1DA3412A12E79CCC9C81264CFE9AD080997058FB61C4243CC32DBAA156B9
    G = (
        0x79A70B2B70400553AE7C9DF416C792C61128751AC92969240C25A07D728BDC93E21F7787ED6972249DE732F38496CD11698713093E9C04FC,
        0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF80000000000000000000000000000000000000000000000000000001,
    )
    n = 0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF7CCA23E9C44EDB49AED63690216CC2728DC58F552378C292AB5844F3


class E521(TwistedEdwardsInterface):
    p = 0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    a = 1
    d = 0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA4331
    G = (
        0x752CB45C48648B189DF90CB2296B2878A3BFD9F42FC6C818EC8BF3C9C0C6203913F6ECC5CCC72434B1AE949D568FC99C6059D0FB13364838AA302A940A2F19BA6C,
        0x0C,
    )
    n = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD15B6C64746FC85F736B8AF5E7EC53F04FBD8C4569A8F1F4540EA2435F5180D6B


class Curve41417(TwistedEdwardsInterface):
    p = 0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEF
    a = 0x01
    d = 0xE21
    G = (
        0x1A334905141443300218C0631C326E5FCD46369F44C03EC7F57FF35498A4AB4D6D6BA111301A73FAA8537C64C4FD3812F3CBC595,
        0x22,
    )
    n = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEB3CC92414CF706022B36F1C0338AD63CF181B0E71A5E106AF79


class JubJub(TwistedEdwardsInterface):
    p = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
    a = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000000
    d = 0x2A9318E74BFA2B48F5FD9207E6BD7FD4292D7F6D37579D2601065FD6D6343EB1
    G = (
        0x11DAFE5D23E1218086A365B99FBF3D3BE72F6AFD7D1F72623E6B071492D1122B,
        0x1D523CF1DDAB1A1793132E78C866C0C33E26BA5CC220FED7CC3F870E59D292AA,
    )
    n = 0xE7DB4EA6533AFA906673B0101343B00A6682093CCC81082D0970E5ED6F72CB7


class MDC201601(TwistedEdwardsInterface):
    p = 0xF13B68B9D456AFB4532F92FDD7A5FD4F086A9037EF07AF9EC13710405779EC13
    a = 1
    d = 0x571304521965B68A7CDFBFCCFB0CB9625F1270F63F21F041EE9309250300CF89
    G = (
        0xB681886A7F903B83D85B421E03CBCF6350D72ABB8D2713E2232C25BFEE68363B,
        0xCA6734E1B59C0B0359814DCF6563DA421DA8BC3D81A93A3A7E73C355BD2864B5,
    )
    n = 0x3C4EDA2E7515ABED14CBE4BF75E97F534FB38975FAF974BB588552F421B0F7FB


class Numsp256t1(TwistedEdwardsInterface):
    p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF43
    a = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF42
    d = 0x3BEE
    G = (0x0D, 0x7D0AB41E2A1276DBA3D330B39FA046BFBE2A6D63824D303F707F6FB5331CADBA)
    n = 0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBE6AA55AD0A6BC64E5B84E6F1122B4AD


class Numsp384t1(TwistedEdwardsInterface):
    p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEC3
    a = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEC2
    d = 0x5158A
    G = (
        0x08,
        0x749CDABA136CE9B65BD4471794AA619DAA5C7B4C930BFF8EBD798A8AE753C6D72F003860FEBABAD534A4ACF5FA7F5BEE,
    )
    n = 0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFECD7D11ED5A259A25A13A0458E39F4E451D6D71F70426E25


class Numsp512t1(TwistedEdwardsInterface):
    p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7
    a = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC6
    d = 0x9BAA8
    G = (
        0x20,
        0x7D67E841DC4C467B605091D80869212F9CEB124BF726973F9FF048779E1D614E62AE2ECE5057B5DAD96B7A897C1D72799261134638750F4F0CB91027543B1C5E,
    )
    n = 0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA7E50809EFDABBB9A624784F449545F0DCEA5FF0CB800F894E78D1CB0B5F0189


class Id_tc26_gost_3410_2012_256_paramSetA(TwistedEdwardsInterface):
    p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97
    a = 0x01
    d = 0x605F6B7C183FA81578BC39CFAD518132B9DF62897009AF7E522C32D6DC7BFFB
    G = (0x0D, 0x60CA1E32AA475B348488C38FAB07649CE7EF8DBE87F22E81F92B2592DBA300E7)
    n = 0x400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67


class Id_tc26_gost_3410_2012_512_paramSetC(TwistedEdwardsInterface):
    p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7
    a = 0x01
    d = 0x9E4F5D8C017D8D9F13A5CF3CDF5BFE4DAB402D54198E31EBDE28A0621050439CA6B39E0A515C06B304E2CE43E79E369E91A0CFC2BC2A22B4CA302DBB33EE7550
    G = (
        0x12,
        0x469AF79D1FB1F5E16B99592B77A01E2A0FDFB0D01794368D9A56117F7B38669522DD4B650CF789EEBF068C5D139732F0905622C04B2BAAE7600303EE73001A3D,
    )
    n = 0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC98CDBA46506AB004C33A9FF5147502CC8EDA9E7A769A12694623CEF47F023ED


class Test_Curve(TwistedEdwardsInterface):
    p = 13
    a = 1
    d = 2
    G = (9, 4)
    n = 8

    def __init__(self):
        logger.warn(
            "edwards test-curve is for development and educational purposes only"
            " and should not be used in production."
        )
