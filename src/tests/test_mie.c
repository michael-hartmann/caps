#include "sfunc.h"
#include "libcasimir.h"
#include "unittest.h"

#include "test_mie.h"

static float64 _mie_lna_perf(int l, double arg, sign_t *sign_a)
{
    float64 lna, lnb;
    sign_t dummy;

    casimir_t self;
    casimir_init(&self, 1, 2*arg);

    casimir_lnab_perf(&self, 1, l, &lna, &lnb, sign_a, &dummy);

    return lna;
}

static float64 _mie_lnb_perf(int l, double arg, sign_t *sign_b)
{
    float64 lna, lnb;
    sign_t dummy;

    casimir_t self;
    casimir_init(&self, 1, 2*arg);

    casimir_lnab_perf(&self, 1, l, &lna, &lnb, &dummy, sign_b);

    return lnb;
}

int test_mie(void)
{
    sign_t sign;
    unittest_t test;
    unittest_init(&test, "Mie", "Test Mie functions al,bl for various parameters");

    AssertAlmostEqual(&test, _mie_lna_perf(3,3,&sign), 1.69245030620195999527136501278285739542);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, _mie_lnb_perf(3,3,&sign), 1.56690522122601151344771063007338242470);
    AssertEqual(&test, sign, +1);

    AssertAlmostEqual(&test, _mie_lna_perf(5,3,&sign), -3.0833346627751969447764584120255355928);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, _mie_lnb_perf(5,3,&sign), -3.206110089012862);
    AssertEqual(&test, sign, +1);

    AssertAlmostEqual(&test, _mie_lna_perf(6,3,&sign), -5.9784534006614519534263858905222671508);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(6,3,&sign), -6.093433624873396);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(50,1,&sign), -365.79392430155136873254652598497788467);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(50,1,&sign), -365.8137152817732);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(50,100,&sign), 174.312311394717498345493328440091268638);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(50,100,&sign), 174.3104974165916);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(100,2,&sign), -726.30666287313622505782248775666927440);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,2,&sign), -726.31660729529160039068920762222369746);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(100,100,&sign), 104.995538882222840217449054833161119273);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,100,&sign), 104.991994545284307924459751266769389806);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(100,200,&sign), 349.797396142286487697982285728323520693);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,200,&sign), 349.796495444169224247914686109960469914);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(100,300,&sign), 565.945090427442791217282609438111518146);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(100,300,&sign), 565.944771508594240527789448563293872224);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(40,0.01,&sign), -648.64173507113317100791502344721329801);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(40,0.01,&sign), -648.66642768146380908837317082847275023);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(4,0.01,&sign), -52.728523465163079623273628437431163764);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(4,0.01,&sign), -52.951665263244193273611377733941258468);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(100, 1,&sign), -865.64416766125947723847038002565294466);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(100, 1,&sign), -865.65411651438956672257219858405454738);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(200, 1,&sign), -2003.2650296097250827819811723946036985);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(200, 1,&sign), -2003.2700169651353702728300966686213123);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(500, 1,&sign), -5915.3540172801973222879731127707420165);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(500, 1,&sign), -5915.3560152708959234646060460452676614);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(1000,1,&sign), -13210.097886016582816381723859188214382);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(1000,1,&sign), -13210.098885515418147666770491500341304);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(1500,1,&sign), -21027.801495754698520339005510266151145);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(1500,1,&sign), -21027.802162198797680563024091196856145);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(2000,1,&sign), -29185.185215718437245666170730719302820);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(2000,1,&sign), -29185.185715593291537268801628970184306);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(4000,1,&sign), -63907.254638502745089438029419051667513);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(4000,1,&sign), -63907.254888471476868081861216637606434);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(6000,1,&sign), -100722.02877381278699710526450383761738);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(6000,1,&sign), -100722.02894046555937519142066187010733);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(8000,1,&sign), -138895.87737532557472806993280328577996);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(8000,1,&sign), -138895.87750031775994991231816543644142);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(8000,2,&sign), -127804.82915169761553869615467279017278);
    AssertAlmostEqual(&test, _mie_lnb_perf(8000,2,&sign), -127804.82927668979197312446971247136798);

    AssertAlmostEqual(&test, _mie_lna_perf(8000,10,&sign), -102052.20711521595799404201856120974620);
    AssertAlmostEqual(&test, _mie_lnb_perf(8000,10,&sign), -102052.20724020785323176373561783658054);

    AssertAlmostEqual(&test, _mie_lna_perf(8000,100,&sign), -65207.924343109906560378428855687546567);
    AssertAlmostEqual(&test, _mie_lnb_perf(8000,100,&sign), -65207.924468072809047480677399738325223);

    AssertAlmostEqual(&test, _mie_lna_perf(8000,1000,&sign), -28302.510543208722891435876522311264098);
    AssertAlmostEqual(&test, _mie_lnb_perf(8000,1000,&sign), -28302.510665327950709002165702496426180);

    AssertAlmostEqual(&test, _mie_lna_perf(8000,0.1,&sign), -175739.54151019442982872552305287434276);
    AssertAlmostEqual(&test, _mie_lnb_perf(8000,0.1,&sign), -175739.54163518661795041477757735386724);

    AssertAlmostEqual(&test, _mie_lna_perf(8000,0.01,&sign), -212583.20558381086668504093590666144975);
    AssertAlmostEqual(&test, _mie_lnb_perf(8000,0.01,&sign), -212583.20570880305483572865968901322558);

    AssertAlmostEqual(&test, _mie_lna_perf(8000,0.001,&sign), -249426.86965681477878812491123551949706);
    AssertAlmostEqual(&test, _mie_lnb_perf(8000,0.001,&sign), -249426.86978180696693910261971050662029);

    AssertAlmostEqual(&test, _mie_lna_perf(3,3,&sign), 1.692450306201961);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, _mie_lnb_perf(3,3,&sign), 1.56690522122601151344771063007338242470);
    AssertEqual(&test, sign, +1);

    AssertAlmostEqual(&test, _mie_lna_perf(4,3,&sign), -0.50863950281017);
    AssertEqual(&test, sign, 1);
    AssertAlmostEqual(&test, _mie_lnb_perf(4,3,&sign), -0.6364737229043445034503593249415023186);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(70,1,&sign), -557.44938197296952054007696337637816079);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(70,1,&sign), -557.46356232717235072070323153088341934);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(70,2,&sign), -459.69436414673188217327579826124991504);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(70,2,&sign), -459.70853167148581678807565076374055502);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(70,70,&sign), 73.0260264952860166234340855388092441193);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(70,70,&sign), 73.0209577995230343081965431939009180221);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(70,100,&sign), 151.413559054452923417338335312392371628);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(70,100,&sign), 151.410845411293235441480855068620711752);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(7,0.2,&sign), -50.341577269323421928601221896653766867);
    AssertEqual(&test, sign, -1);
    AssertAlmostEqual(&test, _mie_lnb_perf(7,0.2,&sign), -50.474963358901354311599799537450688793);
    AssertEqual(&test, sign, +1);

    AssertAlmostEqual(&test, _mie_lna_perf(20,0.1,&sign), -206.31468726371069221722689886250847077);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(20,0.1,&sign), -206.36347568161690202924645177884481860);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(20,0.01,&sign), -300.72091638627787595270816759792424287);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(20,0.01,&sign), -300.76970653298415620564166959511978878);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(30,0.01,&sign), -471.34450706689552653481304006884802111);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(30,0.01,&sign), -471.37729688442461998847047634057026492);
    AssertEqual(&test, sign, -1);

    AssertAlmostEqual(&test, _mie_lna_perf(30,0.001,&sign), -611.80219935898866909805624025148891633);
    AssertEqual(&test, sign, +1);
    AssertAlmostEqual(&test, _mie_lnb_perf(30,0.001,&sign), -611.83498918175872098732849032138015877);
    AssertEqual(&test, sign, -1);

    return test_results(&test, stderr);
}
