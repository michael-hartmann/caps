from unittest import main, TestCase
import numpy as np
import mpmath as mp
import libcasimir
from math import log,e
from scipy.misc import factorial2
import itertools

# increase accuracy of mpmath
mp.mp.dps = 50

class TestBessel(TestCase):
    """Test Bessel functions"""

    def test_I_explicit(self):
        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=10)

        lnInu = libcasimir.sfunc.lnInu

        assertAlmostEqual(lnInu(0,1e-6), -7.133546631626697)
        assertAlmostEqual(lnInu(0,1e-5), -5.982254085113174)
        assertAlmostEqual(lnInu(0,1e-4), -4.830961536966152)
        assertAlmostEqual(lnInu(0,1e-3), -3.679668825469134)
        assertAlmostEqual(lnInu(0,1e-2), -2.528359779027661)
        assertAlmostEqual(lnInu(0,1e-1), -1.375417787678169)
        assertAlmostEqual(lnInu(0,1e0),  -0.064351991073531)
        assertAlmostEqual(lnInu(0,5e0),   3.276297109617906)
        assertAlmostEqual(lnInu(0,1e1),   7.929768918237150)
        assertAlmostEqual(lnInu(0,1e2),   96.77847637380128)
        assertAlmostEqual(lnInu(0,1e3),   995.6271838273042)
        assertAlmostEqual(lnInu(0,1e4),   9994.475891280807)
        assertAlmostEqual(lnInu(0,1e5),   99993.32459873431)
        assertAlmostEqual(lnInu(0,1e6),   999992.17330618785)

        assertAlmostEqual(lnInu(1,1e-6), -22.04766947825915)
        assertAlmostEqual(lnInu(1,1),    -1.225791352644727)
        assertAlmostEqual(lnInu(1,3),     1.131235470744604)
        assertAlmostEqual(lnInu(1,1e6),   999992.1733051878)

        assertAlmostEqual(lnInu(2,1e-6), -37.47261794865755)
        assertAlmostEqual(lnInu(2,1),    -2.862970265776753)
        assertAlmostEqual(lnInu(2,5),     2.622265862896675)
        assertAlmostEqual(lnInu(2,1e6),   999992.1733031878)

        assertAlmostEqual(lnInu(3,1), -4.824473578629219)

        assertAlmostEqual(lnInu(23,1e-6), -394.1439513814884)
        assertAlmostEqual(lnInu(23,5),    -31.40382021014728)
        assertAlmostEqual(lnInu(23,1e6),   999992.1730301876)

        assertAlmostEqual(lnInu(119,1e-6), -2189.202200199878)
        assertAlmostEqual(lnInu(119,0.5),  -621.0792579289692)
        assertAlmostEqual(lnInu(119,3),    -406.9458492626251)
        assertAlmostEqual(lnInu(119,30),   -129.9524456900199)
        assertAlmostEqual(lnInu(119,300),   272.6929318295042)
        assertAlmostEqual(lnInu(119,1e6),   999992.1661661842)

        assertAlmostEqual(lnInu(702,1e-6),  -14098.666835519577122094)
        assertAlmostEqual(lnInu(702,1e-4),  -10863.534779862939382744)
        assertAlmostEqual(lnInu(702,3),     -3621.4923374733442116413)
        assertAlmostEqual(lnInu(702,1234),   1034.4300403851143436433)
        assertAlmostEqual(lnInu(702,12345),  12319.387046237228462572)

        assertAlmostEqual(lnInu(1000,1e-6), -20431.494498396182845040)
        assertAlmostEqual(lnInu(1000,1e-3), -13520.285341774305033047)
        assertAlmostEqual(lnInu(1000,1e-2), -11216.548956209049395583)
        assertAlmostEqual(lnInu(1000,1e-1), -8912.8125681972135967096)
        assertAlmostEqual(lnInu(1000,0.5),  -7302.5698768967633689718)
        assertAlmostEqual(lnInu(1000,1),    -6609.0759355273959800953)
        assertAlmostEqual(lnInu(1000,2),    -5915.5814325009519062869)
        assertAlmostEqual(lnInu(1000,3),    -5509.9123437129452673289)
        assertAlmostEqual(lnInu(1000,1e3),  527.852986878681152219404)
        assertAlmostEqual(lnInu(1000,1e6),  999991.672805979313094335)
        assertAlmostEqual(lnInu(1000,1e10), 9999999987.56808595182509)
        assertAlmostEqual(lnInu(1000,1e15), 999999999999981.811673268)

        assertAlmostEqual(lnInu(2000,1e-6), -42234.894795130034683047337)
        assertAlmostEqual(lnInu(2000,1e-3), -28415.930359526144472047834)
        assertAlmostEqual(lnInu(2000,1e-2), -23809.608880979190355464161)
        assertAlmostEqual(lnInu(2000,1e-1), -19203.287401208029324997325)
        assertAlmostEqual(lnInu(2000,0.5),  -15983.606827406095017386259)
        assertAlmostEqual(lnInu(2000,1),    -14596.965799016187883100221)
        assertAlmostEqual(lnInu(2000,2),    -13210.324489587114945092216)
        assertAlmostEqual(lnInu(2000,3),    -12399.190916285384009011161)
        assertAlmostEqual(lnInu(2000,1e3),  -656.70014078181532009955088)
        assertAlmostEqual(lnInu(2000,1e6),  999990.172305854646890506260)
        assertAlmostEqual(lnInu(2000,1e10), 9999999987.56793590182508883)
        assertAlmostEqual(lnInu(2000,1e15), 999999999999981.811673267338)

        assertAlmostEqual(lnInu(4000,1e-6), -87227.2969461230788510330120)
        assertAlmostEqual(lnInu(4000,1e-3), -59592.8219525549768334749645)
        assertAlmostEqual(lnInu(4000,1e-2), -50381.3302880261118940037035)
        assertAlmostEqual(lnInu(4000,1e-1), -41169.8386228849139406141579)
        assertAlmostEqual(lnInu(4000,0.5),  -34731.2822391979185359387942)
        assertAlmostEqual(lnInu(4000,1),    -31958.3468965104293184241022)
        assertAlmostEqual(lnInu(4000,2),    -29185.4114132506615649536407)
        assertAlmostEqual(lnInu(4000,3),    -27563.3479358811252206864128)
        assertAlmostEqual(lnInu(4000,1e3),  -4261.87314558365887901397646)
        assertAlmostEqual(lnInu(4000,1e6),  999984.1713128587906335733650)
        assertAlmostEqual(lnInu(4000,1e10), 9999999987.567335801825058838)
        assertAlmostEqual(lnInu(4000,1e15), 999999999999981.8116732613379)

        assertAlmostEqual(lnInu(8000,1e-6), -179983.996995002381895872387)
        assertAlmostEqual(lnInu(8000,1e-3), -124718.500885505762638222493)
        assertAlmostEqual(lnInu(8000,1e-2), -106296.668849003806973205933)
        assertAlmostEqual(lnInu(8000,1e-1), -87874.8368121956271975910626)
        assertAlmostEqual(lnInu(8000,0.5),  -74998.5287862680135848981838)
        assertAlmostEqual(lnInu(8000,1),    -69453.0047447650649015166130)
        assertAlmostEqual(lnInu(8000,2),    -63907.4806329627981981416669)
        assertAlmostEqual(lnInu(8000,3),    -60663.5568793227244075509925)
        assertAlmostEqual(lnInu(8000,1e3),  -14156.3252028322219055184688)
        assertAlmostEqual(lnInu(8000,1e6),  999960.1694608923673604843437)
        assertAlmostEqual(lnInu(8000,1e10), 9999999987.564935601824938988)
        assertAlmostEqual(lnInu(8000,1e15), 999999999999981.8116732373359)


    def test_I_mpmath_small(self):
        lnInu = libcasimir.sfunc.lnInu

        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=10)

        for nu in range(20):
            for x in np.logspace(log(1e-6), log(1e6), 100, base=e):
                assertAlmostEqual(lnInu(nu,x), mp.log(mp.besseli(nu+0.5,x)))


    def test_K_mpmath_small(self):
        lnKnu = libcasimir.sfunc.lnInu

        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=10)

        for nu in range(20):
            for x in np.logspace(log(1e-6), log(1e6), 100, base=e):
                assertAlmostEqual(lnKnu(nu,x), mp.log(mp.besselk(nu+0.5,x)))


    def test_K_explicit(self):
        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=10)

        lnKnu = libcasimir.sfunc.lnKnu

        assertAlmostEqual(lnKnu(0,1e-6), 7.133545631626864)
        assertAlmostEqual(lnKnu(0,1e-5), 5.9822440851298415)
        assertAlmostEqual(lnKnu(0,1e-4), 4.830861538632819)
        assertAlmostEqual(lnKnu(0,1e-3), 3.6786689921357962)
        assertAlmostEqual(lnKnu(0,1e-2), 2.5183764456387734)
        assertAlmostEqual(lnKnu(0,1e-1), 1.2770838991417504)
        assertAlmostEqual(lnKnu(0,1e0), -0.7742086473552725)
        assertAlmostEqual(lnKnu(0,5e0), -5.5789276035723227)
        assertAlmostEqual(lnKnu(0,1e1), -10.925501193852295)
        assertAlmostEqual(lnKnu(0,1e2), -102.07679374034932)
        assertAlmostEqual(lnKnu(0,1e3), -1003.2280862868463)
        assertAlmostEqual(lnKnu(0,1e4), -10004.379378833343)
        assertAlmostEqual(lnKnu(0,1e5), -100005.53067137984)
        assertAlmostEqual(lnKnu(0,1e6), -1.0000066819639263e6)

        assertAlmostEqual(lnKnu(1,1e-6), 20.949057189590638)
        assertAlmostEqual(lnKnu(1,1),   -0.0810614667953271)
        assertAlmostEqual(lnKnu(1,3),   -3.0358327192375464)
        assertAlmostEqual(lnKnu(1,1e6), -1.0000066819629263e6)

        assertAlmostEqual(lnKnu(2,1e-6), 35.8631800362233)
        assertAlmostEqual(lnKnu(2,1),    1.17170150170004)
        assertAlmostEqual(lnKnu(2,5),   -5.03660331274696)
        assertAlmostEqual(lnKnu(2,1e6), -1.0000066819609263389096196908e6)

        assertAlmostEqual(lnKnu(3,1), 2.8367092652889516)

        assertAlmostEqual(lnKnu(4,1e15), -1.00000000000001704359684e15)

        assertAlmostEqual(lnKnu(23,1e-6), 390.29380377977833)
        assertAlmostEqual(lnKnu(23,5),    27.5314997887589672718741222750056)
        assertAlmostEqual(lnKnu(23,1e6), -1.00000668168792647542217765299e6)

        assertAlmostEqual(lnKnu(119,1e-6), 2183.7257366479472175742539693253862993069)
        assertAlmostEqual(lnKnu(119,0.5),  615.6027856231534)
        assertAlmostEqual(lnKnu(119,3),    401.4690706673959)
        assertAlmostEqual(lnKnu(119,30),   124.44542144141829)
        assertAlmostEqual(lnKnu(119,300), -279.16349731660983)

        assertAlmostEqual(lnKnu(702,1e-4), 10856.28698728117152647293)
        assertAlmostEqual(lnKnu(702,3), 3614.24453577321548255948381274)
        assertAlmostEqual(lnKnu(702,1234), -1042.3815655681729711090061175312483747)
        assertAlmostEqual(lnKnu(702,12345), -12329.50281632819683895772331427)

        assertAlmostEqual(lnKnu(1000,1e-6),  20423.89309606159906635627)
        assertAlmostEqual(lnKnu(1000,1e-3),  13512.68393943972082096386)
        assertAlmostEqual(lnKnu(1000,1e-2),  11208.94755387441573291358)
        assertAlmostEqual(lnKnu(1000,1e-1),  8905.211165857634910126568)
        assertAlmostEqual(lnKnu(1000,0.5),   7294.968474437304432718383)
        assertAlmostEqual(lnKnu(1000,1),     6601.474532693311622436143)
        assertAlmostEqual(lnKnu(1000,2),     5907.980028168368669895212)
        assertAlmostEqual(lnKnu(1000,3),     5502.310936882873879713131)
        assertAlmostEqual(lnKnu(1000,1e3),  -535.8007129753599475405978)
        assertAlmostEqual(lnKnu(1000,1e6),  -1.0000061814642183370632e6)
        assertAlmostEqual(lnKnu(1000,1e10), -1.000000001128708406232e10)
        assertAlmostEqual(lnKnu(1000,1e15), -1.000000000000017043596e15)

        assertAlmostEqual(lnKnu(2000,1e-6), 42226.60049552117744801)
        assertAlmostEqual(lnKnu(2000,1e-3), 28407.63605991728711208)
        assertAlmostEqual(lnKnu(2000,1e-2), 23801.31458137032062668)
        assertAlmostEqual(lnKnu(2000,1e-1), 19194.99310159792271442)
        assertAlmostEqual(lnKnu(2000,0.5),  15975.31252776600339467)
        assertAlmostEqual(lnKnu(2000,1),    14588.67149928239310903)
        assertAlmostEqual(lnKnu(2000,2),    13202.03018947850774122)
        assertAlmostEqual(lnKnu(2000,3),    12390.89661555209004650)
        assertAlmostEqual(lnKnu(2000,1e3),  648.2943193660655188020)
        assertAlmostEqual(lnKnu(2000,1e6),  -1000004.68096559416710)
        assertAlmostEqual(lnKnu(2000,1e10), -10000000011.2869340123)
        assertAlmostEqual(lnKnu(2000,1e15), -1000000000000017.04359)

        assertAlmostEqual(lnKnu(4000,1e-6), 87218.30962431022872707)
        assertAlmostEqual(lnKnu(4000,1e-3), 59583.83463074212667827)
        assertAlmostEqual(lnKnu(4000,1e-2), 50372.34296621325864582)
        assertAlmostEqual(lnKnu(4000,1e-1), 41160.85130107175139474)
        assertAlmostEqual(lnKnu(4000,0.5),  34722.29491737725786430)
        assertAlmostEqual(lnKnu(4000,1),    31949.35957466633700452)
        assertAlmostEqual(lnKnu(4000,2),    29176.42409131284269294)
        assertAlmostEqual(lnKnu(4000,3),    27554.36061378709545753)
        assertAlmostEqual(lnKnu(4000,1e3),  4252.855518809915147217)
        assertAlmostEqual(lnKnu(4000,1e6),  -999998.679978599250821)
        assertAlmostEqual(lnKnu(4000,1e10), -10000000011.2863339123)
        assertAlmostEqual(lnKnu(4000,1e15), -1000000000000017.04359)

        assertAlmostEqual(lnKnu(8000,1e-6), 179974.31658850311302)
        assertAlmostEqual(lnKnu(8000,1e-3), 124708.82047900649375)
        assertAlmostEqual(lnKnu(8000,1e-2), 106286.98844250453731)
        assertAlmostEqual(lnKnu(8000,1e-1), 87865.156405696280207)
        assertAlmostEqual(lnKnu(8000,0.5),  74988.848379766791829)
        assertAlmostEqual(lnKnu(8000,1),    69443.324338257984503)
        assertAlmostEqual(lnKnu(8000,2),    63897.800226432283229)
        assertAlmostEqual(lnKnu(8000,3),    60653.876472753151824)
        assertAlmostEqual(lnKnu(8000,1e3),  14146.637045201018258)
        assertAlmostEqual(lnKnu(8000,1e6),  -999974.6781506338673)
        assertAlmostEqual(lnKnu(8000,1e10), -10000000011.28393371)
        assertAlmostEqual(lnKnu(8000,1e15), -1000000000000017.043)

class TestDoublefact(TestCase):
    """Test double factorial n!!"""

    def test_explicit(self):
        ln_doublefact = libcasimir.sfunc.ln_doublefact

        self.assertEqual(ln_doublefact(0),   0);
        self.assertEqual(ln_doublefact(1),   0);

        self.assertAlmostEqual(ln_doublefact(2),   0.693147180559945309)
        self.assertAlmostEqual(ln_doublefact(3),   1.098612288668109691)
        self.assertAlmostEqual(ln_doublefact(4),   2.079441541679835928)
        self.assertAlmostEqual(ln_doublefact(5),   2.708050201102210065)

        self.assertAlmostEqual(ln_doublefact(51),  77.07730784751820516)
        self.assertAlmostEqual(ln_doublefact(52),  79.28352845556058002)

        self.assertAlmostEqual(ln_doublefact(100), 183.1351259797702975)
        self.assertAlmostEqual(ln_doublefact(101), 185.2193700926344520)

        self.assertAlmostEqual(ln_doublefact(333), 803.8068699169127948)
        self.assertAlmostEqual(ln_doublefact(334), 806.9389802679216196)

        self.assertAlmostEqual(ln_doublefact(499), 1303.998431529316796)
        self.assertAlmostEqual(ln_doublefact(500), 1307.332026930839287)


    def test_scipy(self):
        ln_doublefact = libcasimir.sfunc.ln_doublefact

        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=10)

        every = 200     # test all integers up to 200
        skip  = 13      # then only test every 13-th integer
        maximum = 15000 # up to 15000

        for i in itertools.chain(range(0, every), range(every, maximum+1, skip)):
            x = ln_doublefact(i)
            y = mp.log(mp.mpf(factorial2(i, exact=True)))
            assertAlmostEqual(x,y)


class TestLambda(TestCase):
    casimir = libcasimir.Casimir()

    def lnLambda_mpmath(self,l1,l2,m):
        num   = (2*l1+1)*(2*l2+1)*mp.factorial(l1-m)*mp.factorial(l2-m)
        denom = mp.factorial(l1+m)*mp.factorial(l2+m)*l1*(l1+1)*l2*(l2+1)
        return mp.log(mp.sqrt( num/denom ))

        
    def test_explicit(self):
        lnLambda = self.casimir.lnLambda

        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=10)

        assertAlmostEqual(lnLambda(1,1,0),           0.40546510810816438197)
        assertAlmostEqual(lnLambda(1,1,1),          -0.28768207245178092743)

        assertAlmostEqual(lnLambda(2,1,1),          -1.13088154923689527723)
        assertAlmostEqual(lnLambda(4,5,3),          -10.1192134441661378305)
        assertAlmostEqual(lnLambda(5,7,3),          -12.0792090317047992966)
        assertAlmostEqual(lnLambda(16,6,4),         -20.0208410616292587591)

        assertAlmostEqual(lnLambda(50,50,0),        -3.22872812131121237937)
        assertAlmostEqual(lnLambda(50,50,1),        -11.0725767594636842096)
        assertAlmostEqual(lnLambda(50,50,50),       -366.968103676874702523)
        assertAlmostEqual(lnLambda(50,20,10),       -71.7083841252765817069)

        assertAlmostEqual(lnLambda(100,1,0),        -1.75576034333105534293)
        assertAlmostEqual(lnLambda(100,1,1),        -6.71247928502570340710)
        assertAlmostEqual(lnLambda(100,6,4),        -28.1900765794504258135)
        assertAlmostEqual(lnLambda(100,33,17),      -140.561511466323128452)
        assertAlmostEqual(lnLambda(100,100,0),      -3.91698579477027506785)
        assertAlmostEqual(lnLambda(100,100,50),     -460.459324692420926858)
        assertAlmostEqual(lnLambda(100,201,10),     -103.403675574453234983)

        assertAlmostEqual(lnLambda(200,200,0),      -4.60766084730054324266)
        assertAlmostEqual(lnLambda(200,100,70),     -689.799517140546177067)
        assertAlmostEqual(lnLambda(500,500,0),      -5.52245942019183595607)

        assertAlmostEqual(lnLambda(1000,1000,0),    -6.21510772371362422788)
        assertAlmostEqual(lnLambda(1000,1000,1),    -20.03161778201098186516)
        assertAlmostEqual(lnLambda(1000,1000,2),    -33.8481258423043454917)
        assertAlmostEqual(lnLambda(1000,1000,3),    -47.6646299065777449973)
        assertAlmostEqual(lnLambda(1000,1000,10),   -144.379878626542707199)
        assertAlmostEqual(lnLambda(1000,1000,20),   -282.542651228911456181)
        assertAlmostEqual(lnLambda(1000,1000,50),   -696.998971043064390454)
        assertAlmostEqual(lnLambda(1000,1000,100),  -1387.53214390965801576)
        assertAlmostEqual(lnLambda(1000,1000,499),  -6855.76100769313908061)
        assertAlmostEqual(lnLambda(1000,1000,500),  -6869.29083418131424684)
        assertAlmostEqual(lnLambda(1000,1000,501),  -6882.81932911136990054)
        assertAlmostEqual(lnLambda(1000,1000,999),  -13205.1385557779782981)
        assertAlmostEqual(lnLambda(1000,1000,1000), -13212.7394582375203805)

        assertAlmostEqual(lnLambda(1500,1500,0),    -6.62040637328339627357)
        assertAlmostEqual(lnLambda(1500,1500,1),    -21.2475135920071596617)
        assertAlmostEqual(lnLambda(1500,1500,2),    -35.8746199224338374194)
        assertAlmostEqual(lnLambda(1500,1500,500),  -7301.01657839892680370)
        assertAlmostEqual(lnLambda(1500,1500,1000), -14460.2359699258235793)
        assertAlmostEqual(lnLambda(1500,1500,1500), -21030.6452594188310976)

        assertAlmostEqual(lnLambda(1500,1,0),    -3.1074706325876159457)

        assertAlmostEqual(lnLambda(5000,1,1),    -12.573207215572773561)
        assertAlmostEqual(lnLambda(5000,500,1),  -21.406202989196329346)
        assertAlmostEqual(lnLambda(5000,2000,1), -23.484521169044081354)
        assertAlmostEqual(lnLambda(5000,5000,1), -24.858732358693766195)


    def test_mpmath_small(self):
        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=10)

        lnLambda = self.casimir.lnLambda
        lnLambda_mpmath = self.lnLambda_mpmath

        for l1 in range(1,10):
            for l2 in range(1,10):
                lmin = min(l1,l2)
                for m in range(lmin+1):
                    assertAlmostEqual(lnLambda(l1,l2,m), lnLambda_mpmath(l1,l2,m))


    def test_mpmath_arbitrary(self):
        casimir = self.casimir
        lmin = 1      # start from l1,l2 = 1
        lmax = 15000 # up to l1,l2=lmax
        N = 50

        assertAlmostEqual = lambda x,y: self.assertAlmostEqual(x,y, places=10)

        lnLambda = self.casimir.lnLambda
        lnLambda_mpmath = self.lnLambda_mpmath

        for l1 in map(int, np.logspace(log(lmin), log(lmax), N, base=e)):
            for l2 in map(int, np.logspace(log(lmin), log(lmax), N, base=e)):
                min_l1l2 = min(l1,l2)
                for m in map(int, np.logspace(log(lmin), log(min_l1l2), N, base=e)):
                    assertAlmostEqual(lnLambda(l1,l2,m), lnLambda_mpmath(l1,l2,m))


if __name__ == "__main__":
    main()
