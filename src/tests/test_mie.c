#include "misc.h"
#include "libcasimir.h"
#include "unittest.h"

#include "test_mie.h"

int test_mie(void)
{
    sign_t sign_a, sign_b;
    double lna, lnb;
    unittest_t test;
    unittest_init(&test, "Mie", "Test Mie functions al,bl for various parameters", 1e-10);

    casimir_t *self = casimir_init(1);

    casimir_lnab_perf(self, 6, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 1.69245030620195999527136501278285739542);
    AssertAlmostEqual(&test, lnb, 1.56690522122601151344771063007338242470);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab_perf(self, 6, 5, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -3.0833346627751969447764584120255355928);
    AssertAlmostEqual(&test, lnb, -3.206110089012862);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab_perf(self, 6, 6, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -5.9784534006614519534263858905222671508);
    AssertAlmostEqual(&test, lnb, -6.093433624873396);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 50, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -365.79392430155136873254652598497788467);
    AssertAlmostEqual(&test, lnb, -365.8137152817732);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 200, 50, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 174.312311394717498345493328440091268638);
    AssertAlmostEqual(&test, lnb, 174.3104974165916);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 4, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -726.30666287313622505782248775666927440);
    AssertAlmostEqual(&test, lnb, -726.31660729529160039068920762222369746);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 200, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 104.995538882222840217449054833161119273);
    AssertAlmostEqual(&test, lnb, 104.991994545284307924459751266769389806);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 400, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 349.797396142286487697982285728323520693);
    AssertAlmostEqual(&test, lnb, 349.796495444169224247914686109960469914);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 600, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 565.945090427442791217282609438111518146);
    AssertAlmostEqual(&test, lnb, 565.944771508594240527789448563293872224);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 0.02, 40, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -648.64173507113317100791502344721329801);
    AssertAlmostEqual(&test, lnb, -648.66642768146380908837317082847275023);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 0.02, 4, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -52.728523465163079623273628437431163764);
    AssertAlmostEqual(&test, lnb, -52.951665263244193273611377733941258468);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 100, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -865.64416766125947723847038002565294466);
    AssertAlmostEqual(&test, lnb, -865.65411651438956672257219858405454738);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 200, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -2003.2650296097250827819811723946036985);
    AssertAlmostEqual(&test, lnb, -2003.2700169651353702728300966686213123);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 500, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -5915.3540172801973222879731127707420165);
    AssertAlmostEqual(&test, lnb, -5915.3560152708959234646060460452676614);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 1000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -13210.097886016582816381723859188214382);
    AssertAlmostEqual(&test, lnb, -13210.098885515418147666770491500341304);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 1500, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -21027.801495754698520339005510266151145);
    AssertAlmostEqual(&test, lnb, -21027.802162198797680563024091196856145);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 2000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -29185.185215718437245666170730719302820);
    AssertAlmostEqual(&test, lnb, -29185.185715593291537268801628970184306);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 4000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -63907.254638502745089438029419051667513);
    AssertAlmostEqual(&test, lnb, -63907.254888471476868081861216637606434);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 6000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -100722.02877381278699710526450383761738);
    AssertAlmostEqual(&test, lnb, -100722.02894046555937519142066187010733);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 8000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -138895.87737532557472806993280328577996);
    AssertAlmostEqual(&test, lnb, -138895.87750031775994991231816543644142);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 4, 8000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -127804.82915169761553869615467279017278);
    AssertAlmostEqual(&test, lnb, -127804.82927668979197312446971247136798);

    casimir_lnab_perf(self, 20, 8000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -102052.20711521595799404201856120974620);
    AssertAlmostEqual(&test, lnb, -102052.20724020785323176373561783658054);

    casimir_lnab_perf(self, 200, 8000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -65207.924343109906560378428855687546567);
    AssertAlmostEqual(&test, lnb, -65207.924468072809047480677399738325223);

    casimir_lnab_perf(self, 2000, 8000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -28302.510543208722891435876522311264098);
    AssertAlmostEqual(&test, lnb, -28302.510665327950709002165702496426180);

    casimir_lnab_perf(self, 0.2, 8000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -175739.54151019442982872552305287434276);
    AssertAlmostEqual(&test, lnb, -175739.54163518661795041477757735386724);

    casimir_lnab_perf(self, 0.02, 8000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -212583.20558381086668504093590666144975);
    AssertAlmostEqual(&test, lnb, -212583.20570880305483572865968901322558);

    casimir_lnab_perf(self, 0.002, 8000, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -249426.86965681477878812491123551949706);
    AssertAlmostEqual(&test, lnb, -249426.86978180696693910261971050662029);

    casimir_lnab_perf(self, 6, 3, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 1.692450306201961);
    AssertAlmostEqual(&test, lnb, 1.56690522122601151344771063007338242470);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab_perf(self, 6, 4, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -0.50863950281017);
    AssertAlmostEqual(&test, lnb, -0.6364737229043445034503593249415023186);
    AssertEqual(&test, sign_a, 1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 2, 70, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -557.44938197296952054007696337637816079);
    AssertAlmostEqual(&test, lnb, -557.46356232717235072070323153088341934);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 4, 70, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -459.69436414673188217327579826124991504);
    AssertAlmostEqual(&test, lnb, -459.70853167148581678807565076374055502);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 140, 70, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 73.0260264952860166234340855388092441193);
    AssertAlmostEqual(&test, lnb, 73.0209577995230343081965431939009180221);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 200, 70, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, 151.413559054452923417338335312392371628);
    AssertAlmostEqual(&test, lnb, 151.410845411293235441480855068620711752);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 0.4, 7, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -50.341577269323421928601221896653766867);
    AssertAlmostEqual(&test, lnb, -50.474963358901354311599799537450688793);
    AssertEqual(&test, sign_a, -1);
    AssertEqual(&test, sign_b, +1);

    casimir_lnab_perf(self, 0.2, 20, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -206.31468726371069221722689886250847077);
    AssertAlmostEqual(&test, lnb, -206.36347568161690202924645177884481860);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 0.02, 20, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -300.72091638627787595270816759792424287);
    AssertAlmostEqual(&test, lnb, -300.76970653298415620564166959511978878);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 0.02, 30, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -471.34450706689552653481304006884802111);
    AssertAlmostEqual(&test, lnb, -471.37729688442461998847047634057026492);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_lnab_perf(self, 0.002, 30, &lna, &lnb, &sign_a, &sign_b);
    AssertAlmostEqual(&test, lna, -611.80219935898866909805624025148891633);
    AssertAlmostEqual(&test, lnb, -611.83498918175872098732849032138015877);
    AssertEqual(&test, sign_a, +1);
    AssertEqual(&test, sign_b, -1);

    casimir_free(self);

    return test_results(&test, stderr);
}
