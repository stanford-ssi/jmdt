#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

/* Logarithm of interpolated densities, from 0 to 1000km.
 * Generated from some old legit-looking FORTRAN script. */
const double us1976[] = {0.20294084399669038,
			 0.10589037525743217,
			 0.0065783153601225068,
			 -0.095135195115038637,
			 -0.19924393598990395,
			 -0.30594109172406203,
			 -0.41534879118234508,
			 -0.52759884435182181,
			 -0.64285338550547833,
			 -0.7612975499086565,
			 -0.88307358101395828,
			 -1.0084060207819607,
			 -1.1649444173584578,
			 -1.3220058712375287,
			 -1.4790238737121673,
			 -1.6359872476676369,
			 -1.7929401659762174,
			 -1.9498177738863303,
			 -2.1066072097526054,
			 -2.2633643798407643,
			 -2.4201306568490892,
			 -2.5807789875735181,
			 -2.740935028437987,
			 -2.9003130087905555,
			 -3.0589276970511383,
			 -3.2167780267860544,
			 -3.3738643555945163,
			 -3.5302360246745588,
			 -3.6858440655704148,
			 -3.8407261235461121,
			 -3.9948612837371504,
			 -4.148251796291011,
			 -4.3009997960315083,
			 -4.4590804801161967,
			 -4.6164940597195345,
			 -4.7720042948477088,
			 -4.9256647481864224,
			 -5.0774965106034724,
			 -5.2275607180972363,
			 -5.3759114079865906,
			 -5.5225364960891801,
			 -5.6675276939970303,
			 -5.8108778276077437,
			 -5.9526670004000701,
			 -6.0928773098510449,
			 -6.231601674798763,
			 -6.3688671174512628,
			 -6.5046262306715201,
			 -6.632626672601976,
			 -6.7569243892795532,
			 -6.8812107237600246,
			 -7.0054783675107517,
			 -7.1238983906516022,
			 -7.2391663449622055,
			 -7.3556061035867391,
			 -7.4732130984011933,
			 -7.5920418311098858,
			 -7.7121022126343082,
			 -7.8334399964981616,
			 -7.9560357353192446,
			 -8.0799710516948817,
			 -8.2052698277164371,
			 -8.3319176685459233,
			 -8.4599625579619655,
			 -8.5894951031716431,
			 -8.7204728428016693,
			 -8.8529456871131007,
			 -8.9869568494573659,
			 -9.122604257672096,
			 -9.2598457918800214,
			 -9.3987323163474183,
			 -9.5393306682496402,
			 -9.6823780685349128,
			 -9.8297910936137907,
			 -9.9785344268590386,
			 -10.128608056734558,
			 -10.280044555497085,
			 -10.432884561887555,
			 -10.587120140727487,
			 -10.742817243274155,
			 -10.900012677120889,
			 -11.058670192692633,
			 -11.218913468823578,
			 -11.380669885849757,
			 -11.544013736253467,
			 -11.709001177996846,
			 -11.875647225539561,
			 -12.04906889672051,
			 -12.224419932348233,
			 -12.401479441872784,
			 -12.580039086578967,
			 -12.759910159889456,
			 -12.940875501358951,
			 -13.122763457425668,
			 -13.305325139085706,
			 -13.488151145859732,
			 -13.670977407231975,
			 -13.853440903999408,
			 -14.035085484401458,
			 -14.215570301413887,
			 -14.394543647443381,
			 -14.571767946853015,
			 -14.74757504395447,
			 -14.922417149457203,
			 -15.096796465392881,
			 -15.271154997218957,
			 -15.445967837622449,
			 -15.621754421393304,
			 -15.798624256017131,
			 -15.97555496417621,
			 -16.151300893389422,
			 -16.324554174742257,
			 -16.493997102467791,
			 -16.658350085415123,
			 -16.816315676667557,
			 -16.966867931101653,
			 -17.110269005670702,
			 -17.247084676285144,
			 -17.377792604056616,
			 -17.503030936040581,
			 -17.623218594127781,
			 -17.7389990415017,
			 -17.850831779736914,
			 -17.958961334294326,
			 -18.063495930019432,
			 -18.16464414933613,
			 -18.262566500858675,
			 -18.357424294317351,
			 -18.449419782292171,
			 -18.538665049862249,
			 -18.625321532035599,
			 -18.709550187717209,
			 -18.791483589798528,
			 -18.871258287372999,
			 -18.949042565307774,
			 -19.024956611609678,
			 -19.099136375769319,
			 -19.171668924202162,
			 -19.242639046591808,
			 -19.312083760318743,
			 -19.380061665758223,
			 -19.446702639979684,
			 -19.512009621174215,
			 -19.576053878442238,
			 -19.638921520076007,
			 -19.70063506930093,
			 -19.761237213398402,
			 -19.820749538494095,
			 -19.879201771069607,
			 -19.936719392127916,
			 -19.993208302366369,
			 -20.048842865384746,
			 -20.103603745648943,
			 -20.157440794961815,
			 -20.210502092550161,
			 -20.262807620544208,
			 -20.314340057209265,
			 -20.365102001956259,
			 -20.415193074846897,
			 -20.464600713565122,
			 -20.513329113119482,
			 -20.561487181383786,
			 -20.60904469285639,
			 -20.656074587406088,
			 -20.702483297763884,
			 -20.748317013813271,
			 -20.793634654613292,
			 -20.83841858484454,
			 -20.882683697380077,
			 -20.926459712953001,
			 -20.969741988587067,
			 -21.012549190765167,
			 -21.054885864596706,
			 -21.096781798081171,
			 -21.138235975068881,
			 -21.179256113473794,
			 -21.219859878465677,
			 -21.260059052409982,
			 -21.299875201258406,
			 -21.339303839194933,
			 -21.378346816699899,
			 -21.417033209785835,
			 -21.455382562205074,
			 -21.493358858758004,
			 -21.531016114531116,
			 -21.568352809814233,
			 -21.605348397490829,
			 -21.642034636696106,
			 -21.67842734320956,
			 -21.714522601310481,
			 -21.750292617915232,
			 -21.785796951043572,
			 -21.821048469874043,
			 -21.856003840343895,
			 -21.890681729299466,
			 -21.925107580263393,
			 -21.959245437407432,
			 -21.993164791388516,
			 -22.026838006495101,
			 -22.060236842735211,
			 -22.093411296326931,
			 -22.126380875796645,
			 -22.159086633216702,
			 -22.191594888420632,
			 -22.22384934446843,
			 -22.25592791057846,
			 -22.287730383737781,
			 -22.319392735320839,
			 -22.35081680480744,
			 -22.381999041105388,
			 -22.41304650206305,
			 -22.443858798746675,
			 -22.474497754998399,
			 -22.504975970321464,
			 -22.535187197998418,
			 -22.565329616084213,
			 -22.59523802829964,
			 -22.624929461882164,
			 -22.654494291465742,
			 -22.683893763088673,
			 -22.713159221984647,
			 -22.742252182926311,
			 -22.771208711566874,
			 -22.799911291324246,
			 -22.828558662216679,
			 -22.857036856505005,
			 -22.8853935772452,
			 -22.913594247594791,
			 -22.941601697495557,
			 -22.969470596504348,
			 -22.997263473088545,
			 -23.024851429607374,
			 -23.052318116155629,
			 -23.079652630908811,
			 -23.106844089149988,
			 -23.133915097172597,
			 -23.160846250727975,
			 -23.187652469764924,
			 -23.214327389346366,
			 -23.240865114281018,
			 -23.267285751226918,
			 -23.293586594583299,
			 -23.319765718811887,
			 -23.345822061729304,
			 -23.37175551159266,
			 -23.397581500318992,
			 -23.423288346361392,
			 -23.448879374429211,
			 -23.474359247641587,
			 -23.499718012720916,
			 -23.52499504893261,
			 -23.550150250998318,
			 -23.575192182157046,
			 -23.600131217984686,
			 -23.6249796737279,
			 -23.649714614582155,
			 -23.674368972357197,
			 -23.69890140880128,
			 -23.723367643018346,
			 -23.747726857333273,
			 -23.77197703529129,
			 -23.79615997179955,
			 -23.820233811581279,
			 -23.844243923538151,
			 -23.868147981906287,
			 -23.891971497902038,
			 -23.915718839583509,
			 -23.939370505018825,
			 -23.962931276022168,
			 -23.986407078112627,
			 -24.009831810784561,
			 -24.033133678320368,
			 -24.056374755976037,
			 -24.079537655255574,
			 -24.102633236021919,
			 -24.125643915356729,
			 -24.148551014137656,
			 -24.171429047889312,
			 -24.194199223293463,
			 -24.216907508999846,
			 -24.239537426454032,
			 -24.262105908187117,
			 -24.284596760113683,
			 -24.30699279703752,
			 -24.329349453651947,
			 -24.351652442631526,
			 -24.373848081128003,
			 -24.395996389320977,
			 -24.418082883727102,
			 -24.440092193111298,
			 -24.462050076188987,
			 -24.483942089880035,
			 -24.505752921025522,
			 -24.527511250201126,
			 -24.549248770282173,
			 -24.570859702218868,
			 -24.592468845191629,
			 -24.614016707099019,
			 -24.635488862377223,
			 -24.656870050408163,
			 -24.678248536814692,
			 -24.699560732394698,
			 -24.7207922659885,
			 -24.742039212321437,
			 -24.763179033680473,
			 -24.784252995514514,
			 -24.805306404712347,
			 -24.826328289788236,
			 -24.847245204298204,
			 -24.868104013385498,
			 -24.888956036239062,
			 -24.909791479714222,
			 -24.930532767033629,
			 -24.951164649197121,
			 -24.971811080245811,
			 -24.992392449604843,
			 -25.012968564474786,
			 -25.03345564800161,
			 -25.053915180760544,
			 -25.074336587499847,
			 -25.094629523452038,
			 -25.114938848672097,
			 -25.135256956139767,
			 -25.15540751305214,
			 -25.175629102690745,
			 -25.195742481209038,
			 -25.215821878565368,
			 -25.235856841641411,
			 -25.255836330520157,
			 -25.275748700655495,
			 -25.295678458696329,
			 -25.315519797667957,
			 -25.335319660900829,
			 -25.355088053931929,
			 -25.374815096121296,
			 -25.394501041910619,
			 -25.41413553352454,
			 -25.433740984579927,
			 -25.453296545016158,
			 -25.472814605613056,
			 -25.492297029084988,
			 -25.511734056326588,
			 -25.531139895011258,
			 -25.550492665148475,
			 -25.5698199225385,
			 -25.589099734162634,
			 -25.608348358121923,
			 -25.627556599687129,
			 -25.64672851705242,
			 -25.665854707296113,
			 -25.684953836361739,
			 -25.704002580696805,
			 -25.723035377732362,
			 -25.742013951153218,
			 -25.760974324969599,
			 -25.779892714893659,
			 -25.798759703391163,
			 -25.817614316686328,
			 -25.836432319690335,
			 -25.855204731621598,
			 -25.873956616497381,
			 -25.892662423427609,
			 -25.911348585595395,
			 -25.92998901033345,
			 -25.948611591038624,
			 -25.967189719664116,
			 -25.985752846493622,
			 -26.004273877971887,
			 -26.022763900242609,
			 -26.041214645443674,
			 -26.059658972346678,
			 -26.078047740475576,
			 -26.096414514117452,
			 -26.114773631560794,
			 -26.13307357952781,
			 -26.151373033075448,
			 -26.169619282321854,
			 -26.187850325166959,
			 -26.206059079167819,
			 -26.224238088381341,
			 -26.242404453439089,
			 -26.260525906717632,
			 -26.278619685831842,
			 -26.29670449390137,
			 -26.314746939515466,
			 -26.332765932531824,
			 -26.350754159454187,
			 -26.368732234779777,
			 -26.386693626254281,
			 -26.404631457318612,
			 -26.422538496124783,
			 -26.440407144805601,
			 -26.458260379735378,
			 -26.476123016249847,
			 -26.493925565444744,
			 -26.511724236067607,
			 -26.529513024148415,
			 -26.547251774130213,
			 -26.565000891861096,
			 -26.582720147801808,
			 -26.600402198414308,
			 -26.618075654017531,
			 -26.635734377827895,
			 -26.653371913923454,
			 -26.670981477124379,
			 -26.688594910449766,
			 -26.706167153393906,
			 -26.723730766251691,
			 -26.741279719287931,
			 -26.758807671295322,
			 -26.776307959787214,
			 -26.793816882375179,
			 -26.81128533826903,
			 -26.82870567836795,
			 -26.846161127445402,
			 -26.863555271126305,
			 -26.88097388823579,
			 -26.898364802273889,
			 -26.915721130727469,
			 -26.933035658728354,
			 -26.950351456879716,
			 -26.967663264179354,
			 -26.984965542902454,
			 -27.00225246930481,
			 -27.019463672408559,
			 -27.036700289675498,
			 -27.053902261428725,
			 -27.071119654515964,
			 -27.088289931151575,
			 -27.105464527562564,
			 -27.122638473813698,
			 -27.139745351407274,
			 -27.156838729861395,
			 -27.173912774668544,
			 -27.191025764923562,
			 -27.208043582735751,
			 -27.225089450008742,
			 -27.242159756987164,
			 -27.259112791349228,
			 -27.276147793936445,
			 -27.29312132731615,
			 -27.310025164446834,
			 -27.326998283888717,
			 -27.34388919478225,
			 -27.360841318612909,
			 -27.377697596113205,
			 -27.394606022421097,
			 -27.411403774390898,
			 -27.42824355056721,
			 -27.445122630644249,
			 -27.461869233183236,
			 -27.478643435507667,
			 -27.495441821686725,
			 -27.512171982090166,
			 -27.528915892991396,
			 -27.545669448829045,
			 -27.562334958140102,
			 -27.579093023068189,
			 -27.595750501846627,
			 -27.612395651805382,
			 -27.629023113265877,
			 -27.64564756279071,
			 -27.662253818359648,
			 -27.678846745411153,
			 -27.695431645257422,
			 -27.711992588103673,
			 -27.728534721493425,
			 -27.745063655969226,
			 -27.761574091301732,
			 -27.778072058640582,
			 -27.794564099002013,
			 -27.811033341282155,
			 -27.827474051989714,
			 -27.843917350791209,
			 -27.860333656388288,
			 -27.876742642256005,
			 -27.893139410649077,
			 -27.909518827374836,
			 -27.925875514796758,
			 -27.942217495445398,
			 -27.958553432980288,
			 -27.974864452480542,
			 -27.991173527827794,
			 -28.007446974214883,
			 -28.023723021336355,
			 -28.039982328293423,
			 -28.056219653727101,
			 -28.072445062890981,
			 -28.08866938792492,
			 -28.104856075355325,
			 -28.121047876779752,
			 -28.137207517546113,
			 -28.153362795655788,
			 -28.169509174511283,
			 -28.185641899694211,
			 -28.201755992448447,
			 -28.217846243256123,
			 -28.233925479562668,
			 -28.250007468588805,
			 -28.266050343153847,
			 -28.282105472116108,
			 -28.298130247939721,
			 -28.314138760851314,
			 -28.330146128175031,
			 -28.34612764784292,
			 -28.362098642523627,
			 -28.378054423046297,
			 -28.393990082859055,
			 -28.409922282507402,
			 -28.425846710521526,
			 -28.441736355254065,
			 -28.457631120035359,
			 -28.473503920657983,
			 -28.489349590464858,
			 -28.505186702752056,
			 -28.521035124724484,
			 -28.536841850019588,
			 -28.552625942626474,
			 -28.568407811097437,
			 -28.584183322578607,
			 -28.599921798075002,
			 -28.615670991566066,
			 -28.63140022818499,
			 -28.647104559157068,
			 -28.662778817611418,
			 -28.678474619255045,
			 -28.694131134754162,
			 -28.709771955501104,
			 -28.725392410339403,
			 -28.741017963463296,
			 -28.756614132928433,
			 -28.772175627862332,
			 -28.787760522773304,
			 -28.803301473521994,
			 -28.8188253390929,
			 -28.834360787392068,
			 -28.849870699484914,
			 -28.865350040269384,
			 -28.880793560501409,
			 -28.896266668726298,
			 -28.911694995184757,
			 -28.927109209162289,
			 -28.942504894640262,
			 -28.957915133766992,
			 -28.973260314852268,
			 -28.988650277834672,
			 -29.003964992680302,
			 -29.019277400361979,
			 -29.034583718614957,
			 -29.049879991921838,
			 -29.065162086505659,
			 -29.080383080516015,
			 -29.095623024317199,
			 -29.110835260569143,
			 -29.126059488054153,
			 -29.141247378291702,
			 -29.156393653244979,
			 -29.171539492898052,
			 -29.186681322682176,
			 -29.201815404770443,
			 -29.216888996850216,
			 -29.231945372597444,
			 -29.247030573745686,
			 -29.262091328489763,
			 -29.277071146859381,
			 -29.292120811686541,
			 -29.307079715704543,
			 -29.322048765402602,
			 -29.337024996338979,
			 -29.351949404396013,
			 -29.366872953218309,
			 -29.381734615720603,
			 -29.396645111551976,
			 -29.411484625414751,
			 -29.426306794933161,
			 -29.44116861346988,
			 -29.455944477935038,
			 -29.470753083705699,
			 -29.485528710010552,
			 -29.500266564048449,
			 -29.514961665702316,
			 -29.529675610906381,
			 -29.544405993280463,
			 -29.559081523721066,
			 -29.573696502473933,
			 -29.588386627829713,
			 -29.602936488671375,
			 -29.617555362420369,
			 -29.632169171238477,
			 -29.646699436938565,
			 -29.661215455039773,
			 -29.675713465281564,
			 -29.690267950322866,
			 -29.704719182981318,
			 -29.719140198802073,
			 -29.733608467115697,
			 -29.748039921963752,
			 -29.762429863233638,
			 -29.776773414895807,
			 -29.791152241127662,
			 -29.805476881212066,
			 -29.819831194663749,
			 -29.834122930437644,
			 -29.848346364971771,
			 -29.862681870585895,
			 -29.876847875103984,
			 -29.891121597316438,
			 -29.905310339767748,
			 -29.919407489522779,
			 -29.933606208922594,
			 -29.947725417483479,
			 -29.961841087937376,
			 -29.975918897308393,
			 -29.989986102358799,
			 -30.004028673219139,
			 -30.018053603000212,
			 -30.032046281735845,
			 -30.046024908392326,
			 -30.059974744737428,
			 -30.073903294461317,
			 -30.087806889999893,
			 -30.101693554488726,
			 -30.115547837898262,
			 -30.129389976005793,
			 -30.143204414053443,
			 -30.156987288444395,
			 -30.170759946160171,
			 -30.184506437281847,
			 -30.198223008443378,
			 -30.211918977370427,
			 -30.225590857527756,
			 -30.239235031875747,
			 -30.252861511213482,
			 -30.266453025328765,
			 -30.280019696438199,
			 -30.293572276959541,
			 -30.307093113804445,
			 -30.320593095784986,
			 -30.334068850257967,
			 -30.34751688061429,
			 -30.360948895550077,
			 -30.374330686675897,
			 -30.387704993445038,
			 -30.401053109386641,
			 -30.414371538092375,
			 -30.427656659326882,
			 -30.440921335035462,
			 -30.454162354368226,
			 -30.4673763906837,
			 -30.480577279324937,
			 -30.493727124413461,
			 -30.506857032907146,
			 -30.519963895301021,
			 -30.533044491045189,
			 -30.546095486465084,
			 -30.559132124885064,
			 -30.572113700524216,
			 -30.585093346223246,
			 -30.598029887463422,
			 -30.610939097188183,
			 -30.623817702830447,
			 -30.636682520560424,
			 -30.649510376605157,
			 -30.66229766900743,
			 -30.675082659879823,
			 -30.68782060923099,
			 -30.700529109145926,
			 -30.713204923213656,
			 -30.7258667923799,
			 -30.738489743222619,
			 -30.75109287173817,
			 -30.763650439583362,
			 -30.776181902974468,
			 -30.788707772484326,
			 -30.801178207749889,
			 -30.813613228719504,
			 -30.82603394053222,
			 -30.838413132855464,
			 -30.850772323832704,
			 -30.863083502600926,
			 -30.875394252070276,
			 -30.887650730408666,
			 -30.89987531060299,
			 -30.912091544389867,
			 -30.924243402925406,
			 -30.936381184018526,
			 -30.94850254455087,
			 -30.960549210900343,
			 -30.972601442022235,
			 -30.984600448651449,
			 -30.996600419988958,
			 -31.008541034236348,
			 -31.020448041991788,
			 -31.032318502591036,
			 -31.044179748989979,
			 -31.055999026256107,
			 -31.067773139004569,
			 -31.079530255141954,
			 -31.091236286055576,
			 -31.102919970763686,
			 -31.114546314137957,
			 -31.126177590524087,
			 -31.137745693805631,
			 -31.149314269466952,
			 -31.160813522721444,
			 -31.172308493769627,
			 -31.183727669274511,
			 -31.195172815974111,
			 -31.206536171087507,
			 -31.217885644635754,
			 -31.229182737121612,
			 -31.240461121216057,
			 -31.251718693990114,
			 -31.262915504559828,
			 -31.274086257903193,
			 -31.285228648511207,
			 -31.296340305633471,
			 -31.307418792540428,
			 -31.318461605848231,
			 -31.329466174912184,
			 -31.340470691787193,
			 -31.351391237133857,
			 -31.362307152692978,
			 -31.373174772501979,
			 -31.383991191502094,
			 -31.39479654489649,
			 -31.405589185902848,
			 -31.416323368577537,
			 -31.427040445600387,
			 -31.437693605748869,
			 -31.448325025298864,
			 -31.458932778738042,
			 -31.469514888068577,
			 -31.480069322195305,
			 -31.490546552398939,
			 -31.501038827211666,
			 -31.511448559766944,
			 -31.521869886125288,
			 -31.532203098178819,
			 -31.542544246315042,
			 -31.552841963836805,
			 -31.563093779162397,
			 -31.573348701621551,
			 -31.583553664694904,
			 -31.593706036158924,
			 -31.603856260753453,
			 -31.613949534872468,
			 -31.623983063862045,
			 -31.63400875673079,
			 -31.644025351210807,
			 -31.653975682063976,
			 -31.663913163837684,
			 -31.673836411606718,
			 -31.683686451159161,
			 -31.693518234711529,
			 -31.703271560677582,
			 -31.713061683694484,
			 -31.722769052303384,
			 -31.732450717873878,
			 -31.742105002154947,
			 -31.751730187581977,
			 -31.761324517004915,
			 -31.770886193452458,
			 -31.780413379934885,
			 -31.78990419928822,
			 -31.799356734062417,
			 -31.808769026456432,
			 -31.818139078303005,
			 -31.827531302911957,
			 -31.836811337464422,
			 -31.846110595852224,
			 -31.855360474766837,
			 -31.864558782344925,
			 -31.873772883350842,
			 -31.882932185279078,
			 -31.892034382536881,
			 -31.901077129904998,
			 -31.910130216595743,
			 -31.91912034457107,
			 -31.928118539817124,
			 -31.937050136375966,
			 -31.945987411921351,
			 -31.954854317024701,
			 -31.963724395271338,
			 -31.972520199426363,
			 -31.981316552440976,
			 -31.99003459629343,
			 -31.998750445901607,
			 -32.007383820977992,
			 -32.016092380788656,
			 -32.024634619854353,
			 -32.033250457919948,
			 -32.041776497435301,
			 -32.050292816632343,
			 -32.058714796785885,
			 -32.067208309561863,
			 -32.075604237268692,
			 -32.08398536902731,
			 -32.09226413311,
			 -32.100612007674819,
			 -32.108854053455445,
			 -32.117075821026475,
			 -32.125276241345482,
			 -32.133454224465524,
			 -32.141517688313961,
			 -32.149646700516456,
			 -32.157749875610079,
			 -32.165732838268077,
			 -32.17378004179524,
			 -32.181703116632569,
			 -32.189594016539857,
			 -32.197547678324433,
			 -32.205371159576757,
			 -32.213256330634145,
			 -32.221007166591598,
			 -32.228818547587224,
			 -32.236561370383527,
			 -32.244263797007527,
			 -32.251934579532602,
			 -32.259582761374013,
			 -32.267197050480299,
			 -32.274776217232208,
			 -32.282329485696692,
			 -32.289845275252276,
			 -32.297332932359005,
			 -32.304791393408898,
			 -32.312230368202229,
			 -32.319627260158676,
			 -32.326991666243138,
			 -32.334333488068047,
			 -32.341640710414204,
			 -32.348923365237937,
			 -32.356169244620091,
			 -32.36338851218364,
			 -32.370580204955289,
			 -32.37774334556746,
			 -32.384876942312843,
			 -32.391979989212615,
			 -32.399051466098939,
			 -32.406102190634122,
			 -32.413131364958183,
			 -32.420126159380565,
			 -32.427085499015355,
			 -32.434020482752466,
			 -32.440930261554747,
			 -32.447801617417177,
			 -32.454658305346946,
			 -32.461487167798992,
			 -32.468274694800321,
			 -32.475057813696083,
			 -32.481797765950525,
			 -32.488519100626817,
			 -32.495208175154659,
			 -32.501877076496527,
			 -32.508511967255721,
			 -32.515125085810723,
			 -32.521715703814621,
			 -32.528269690687523,
			 -32.534812994406806,
			 -32.541317983974793,
			 -32.547810938746117,
			 -32.55426385569249,
			 -32.560689514160153,
			 -32.567101037969344,
			 -32.573497886509337,
			 -32.579851308777876,
			 -32.586188591102129,
			 -32.592494886246456,
			 -32.598783719619412,
			 -32.605040046787082,
			 -32.611292107043383,
			 -32.61749564615797,
			 -32.623693709650631,
			 -32.629856305909051,
			 -32.635997452173129,
			 -32.642131522873044,
			 -32.64822795540784,
			 -32.654301017182888,
			 -32.660350052472054,
			 -32.666389773445481,
			 -32.67238885071059,
			 -32.678377449965943,
			 -32.684339510328535,
			 -32.690274339027567,
			 -32.696197076685998,
			 -32.702107296934024,
			 -32.707972510120136,
			 -32.713823958326863,
			 -32.719661195672778,
			 -32.725483770470213,
			 -32.731274818734022,
			 -32.737033592073729,
			 -32.742775931425378,
			 -32.748518051808901,
			 -32.75422618337835,
			 -32.75991643442223,
			 -32.765588311268068,
			 -32.771224239382896,
			 -32.776857767693649,
			 -32.782471406667035,
			 -32.788064640120155,
			 -32.793636946317854,
			 -32.799187798011538,
			 -32.804716662482008,
			 -32.810223001586422,
			 -32.815724123692419,
			 -32.821201823872009,
			 -32.826655548342849,
			 -32.832084738072695,
			 -32.837507073858973,
			 -32.842903938489798,
			 -32.848293200625463,
			 -32.853656036172339,
			 -32.859010502643557,
			 -32.864337568259174,
			 -32.869655480315409,
			 -32.874944998646257,
			 -32.880224561198332,
			 -32.885474718830196,
			 -32.890714100939817,
			 -32.895942391858156,
			 -32.901159272211274,
			 -32.906344873351166,
			 -32.911518210906671,
			 -32.916678953529832,
			 -32.921806916061236,
			 -32.926941357876863,
			 -32.932042134050903,
			 -32.937148904568559,
			 -32.942221108307599,
			 -32.94727844560925,
			 -32.952320558868507,
			 -32.957347086895901,
			 -32.962357664933819,
			 -32.967351924674801,
			 -32.972329494281681,
			 -32.977310981025681,
			 -32.982254144823202,
			 -32.987200673321247,
			 -32.992107902729664,
			 -32.997017931166724,
			 -33.001909175311418,
			 -33.006802850260321,
			 -33.011655440965804,
			 -33.016509872267584,
			 -33.021344128272609,
			 -33.026157800457298,
			 -33.030972616393598,
			 -33.035743988391033,
			 -33.04053823556881,
			 -33.045310658729008,
			 -33.050060833691646,
			 -33.054788333193997,
			 -33.059515507288531,
			 -33.064219356670478,
			 -33.068922440691864,
			 -33.073601539687424,
			 -33.07827942325828,
			 -33.082932650518593,
			 -33.087584202212952,
			 -33.092210415449223,
			 -33.096834482951976,
			 -33.101456275187751,
			 -33.106051794912986,
			 -33.110644555049213,
			 -33.115210334115559,
			 -33.119772859381087,
			 -33.124331990374451,
			 -33.128863168192161,
			 -33.133390443535276,
			 -33.137913669615401,
			 -33.142407948219322,
			 -33.146897655496893,
			 -33.151382638347492,
			 -33.155862741950273,
			 -33.160312612905983,
			 -33.164782373753638,
			 -33.169221357795095,
			 -33.173654599023592,
			 -33.178056284910227,
			 -33.182477431326198,
			 -33.186866460652624,
			 -33.191274838509948,
			 -33.195650529192172,
			 -33.200019233248454,
			 -33.204380774581928,
			 -33.208734975341592,
			 -33.213055094436527,
			 -33.217393957967936,
			 -33.221724936532958,
			 -33.226047845247884,
			 -33.230335472995655,
			 -33.23464156368577,
			 -33.238939019029388,
			 -33.243200272711753,
			 -33.247479762413072,
			 -33.251750035449462,
			 -33.256010894619649,
			 -33.260234296285027,
			 -33.264475610751717,
			 -33.2687069090099};

double density_us1976(double h);
#endif
