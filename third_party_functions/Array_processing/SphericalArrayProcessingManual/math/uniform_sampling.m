function [a,th,ph]=uniform_sampling(N);
% UNIFORM_SAMPLING generates (nearly) uniform distribution
% of samples on a sphere, having equal sampling weights.
% [a,th,ph]=uniform_sampling(N);
% N is the order.
% a are the sampling weights.
% th are the elevation angles for all samples.
% ph are the azimuth angles for all samples.
%
% Here is the citation for the t-design sampling configurations:
% Hardin and Sloane, McLaren's Improved Snub Cube 
% and Other New Spherical Designs in Three Dimensions.
% Discrete Comput. Geom. 1996,15(4), 429-441. 
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

if N==2 % Total 12 samples
    [a,th,ph]=design3_12_5;

elseif N==3 % Total 32 samples
    [a,th,ph]=design3_32_7;
    
elseif N==4 % Total 36 samples
    [a,th,ph]=design3_36_8;
    
elseif N==6 % Total 84 samples
    [a,th,ph]=design3_84_12;
    
elseif N==8 % Total 144 samples
    [a,th,ph]=design3_144_16;

elseif N==10 % Total 240 samples
    [a,th,ph]=design3_240_21;
    
else
    fprintf('\nNot available for given N');
    return;
end;



function [a,th,ph] = design3_12_5;
% 12 samples, N=2.
% Uses regular isocahedron configuration

j=1:12;

xyz=[
0.850650808352E+00
0			
-0.525731112119E+00
0.525731112119E+00
-0.850650808352E+00
0.000000000000E+00
0		
-0.525731112119E+00
0.850650808352E+00
0.850650808352E+00
0			
0.525731112119E+00
-0.525731112119E+00
-0.850650808352E+00
0
0		
0.525731112119E+00
-0.850650808352E+00
-0.850650808352E+00
0			
-0.525731112119E+00
-0.525731112119E+00
0.850650808352E+00
0
0		
0.525731112119E+00
0.850650808352E+00
-0.850650808352E+00
0			
0.525731112119E+00
0.525731112119E+00
0.850650808352E+00
0
0		
-0.525731112119E+00
-0.850650808352E+00
]';

x=xyz(3*(j-1)+1);
y=xyz(3*(j-1)+2);
z=xyz(3*(j-1)+3);

th=acos(z);
ph=atan2(y,x);

a=ones(size(th))*(4*pi/12);



function [a,th,ph] = design3_32_7;
% 32 samples, N=3.

j=1:32;

xyz=[
.89890353517914318903
.43551221128474708981
.047974454288805545355
.89890353517914318903
-.43551221128474708981
-.047974454288805545355
.89890353517914318903
.047974454288805545355
-.43551221128474708981
.89890353517914318903
-.047974454288805545355
.43551221128474708981
-.89890353517914318903
.43551221128474708981
-.047974454288805545355
-.89890353517914318903
-.43551221128474708981
.047974454288805545355
-.89890353517914318903
.047974454288805545355
.43551221128474708981
-.89890353517914318903
-.047974454288805545355
-.43551221128474708981
.047974454288805545355
.89890353517914318903
.43551221128474708981
-.047974454288805545355
.89890353517914318903
-.43551221128474708981
-.43551221128474708981
.89890353517914318903
.047974454288805545355
.43551221128474708981
.89890353517914318903
-.047974454288805545355
-.047974454288805545355
-.89890353517914318903
.43551221128474708981
.047974454288805545355
-.89890353517914318903
-.43551221128474708981
.43551221128474708981
-.89890353517914318903
.047974454288805545355
-.43551221128474708981
-.89890353517914318903
-.047974454288805545355
.43551221128474708981
.047974454288805545355
.89890353517914318903
-.43551221128474708981
-.047974454288805545355
.89890353517914318903
.047974454288805545355
-.43551221128474708981
.89890353517914318903
-.047974454288805545355
.43551221128474708981
.89890353517914318903
.43551221128474708981
-.047974454288805545355
-.89890353517914318903
-.43551221128474708981
.047974454288805545355
-.89890353517914318903
.047974454288805545355
.43551221128474708981
-.89890353517914318903
-.047974454288805545355
-.43551221128474708981
-.89890353517914318903
.57735026918962576449
.57735026918962576449
.57735026918962576449
.57735026918962576449
.57735026918962576449
-.57735026918962576449
.57735026918962576449
-.57735026918962576449
.57735026918962576449
.57735026918962576449
-.57735026918962576449
-.57735026918962576449
-.57735026918962576449
.57735026918962576449
.57735026918962576449
-.57735026918962576449
.57735026918962576449
-.57735026918962576449
-.57735026918962576449
-.57735026918962576449
.57735026918962576449
-.57735026918962576449
-.57735026918962576449
-.57735026918962576449
]';

x=xyz(3*(j-1)+1);
y=xyz(3*(j-1)+2);
z=xyz(3*(j-1)+3);

th=acos(z);
ph=atan2(y,x);

a=ones(size(th))*(4*pi/32);



function [a,th,ph] = design3_36_8;
% 36 samples, N=4.

j=1:36;

xyz=[
.5074754464108170E+00
-.3062000132395710E+00
.8054254920116630E+00
-.3062000132395690E+00
.8054254920116630E+00
.5074754464108170E+00
-.5074754464108170E+00
.3062000132395700E+00
.8054254920116630E+00
.8054254920116630E+00
.5074754464108170E+00
-.3062000132395690E+00
.3062000132395690E+00
.8054254920116640E+00
-.5074754464108170E+00
.8054254920116630E+00
-.5074754464108170E+00
.3062000132395690E+00
.3062000132395690E+00
-.8054254920116630E+00
.5074754464108160E+00
-.8054254920116630E+00
-.5074754464108170E+00
-.3062000132395690E+00
-.3062000132395700E+00
-.8054254920116640E+00
-.5074754464108160E+00
-.8054254920116630E+00
.5074754464108180E+00
.3062000132395690E+00
.5074754464108170E+00
.3062000132395700E+00
-.8054254920116630E+00
-.5074754464108170E+00
-.3062000132395700E+00
-.8054254920116630E+00
.6263636702652710E+00
-.2435277754091940E+00
-.7405152092807200E+00
-.2435277754091950E+00
-.7405152092807200E+00
.6263636702652710E+00
-.6263636702652710E+00
.2435277754091940E+00
-.7405152092807200E+00
-.7405152092807200E+00
.6263636702652700E+00
-.2435277754091950E+00
.2435277754091950E+00
-.7405152092807190E+00
-.6263636702652710E+00
-.7405152092807200E+00
-.6263636702652700E+00
.2435277754091950E+00
.2435277754091950E+00
.7405152092807190E+00
.6263636702652710E+00
.7405152092807200E+00
-.6263636702652700E+00
-.2435277754091950E+00
-.2435277754091950E+00
.7405152092807200E+00
-.6263636702652710E+00
.7405152092807200E+00
.6263636702652700E+00
.2435277754091950E+00
.6263636702652710E+00
.2435277754091940E+00
.7405152092807200E+00
-.6263636702652710E+00
-.2435277754091940E+00
.7405152092807200E+00
-.2862487234260350E+00
.9571203270924580E+00
-.4452356458542100E-01
.9571203270924580E+00
-.4452356458542000E-01
-.2862487234260350E+00
.2862487234260350E+00
-.9571203270924580E+00
-.4452356458542100E-01
-.4452356458542000E-01
-.2862487234260350E+00
.9571203270924580E+00
-.9571203270924580E+00
-.4452356458541900E-01
.2862487234260350E+00
-.4452356458542100E-01
.2862487234260340E+00
-.9571203270924580E+00
-.9571203270924580E+00
.4452356458542000E-01
-.2862487234260340E+00
.4452356458542100E-01
.2862487234260340E+00
.9571203270924580E+00
.9571203270924580E+00
.4452356458542000E-01
.2862487234260340E+00
.4452356458542100E-01
-.2862487234260340E+00
-.9571203270924580E+00
-.2862487234260340E+00
-.9571203270924580E+00
.4452356458542100E-01
.2862487234260350E+00
.9571203270924580E+00
.4452356458542100E-01
]';

x=xyz(3*(j-1)+1);
y=xyz(3*(j-1)+2);
z=xyz(3*(j-1)+3);

th=acos(z);
ph=atan2(y,x);

a=ones(size(th))*(4*pi/36);



function [a,th,ph] = design3_84_12;
% 84 samples, N=6.

j=1:84;

xyz=[
-.8938049777611360E+00
-.4268621911244970E+00
.1374821134468340E+00
-.4268621912410920E+00
.1374821134452880E+00
-.8938049777056910E+00
.8938049777701570E+00
.4268621911281570E+00
.1374821133768230E+00
.1374821132964000E+00
-.8938049777394910E+00
-.4268621912182710E+00
.4268621912727310E+00
.1374821133773450E+00
.8938049777010320E+00
.1374821135290330E+00
.8938049777077750E+00
.4268621912097560E+00
.4268621911859830E+00
-.1374821134749930E+00
-.8938049777274411E+00
-.1374821133242910E+00
.8938049777252790E+00
-.4268621912390470E+00
-.4268621912174140E+00
-.1374821133472880E+00
.8938049777320730E+00
-.1374821135010710E+00
-.8938049776936550E+00
.4268621912483280E+00
-.8938049776725480E+00
.4268621913280710E+00
-.1374821133907030E+00
.8938049776635530E+00
-.4268621913332600E+00
-.1374821134330650E+00
.9830866003855740E+00
.2230038010752200E-01
-.1817785168533230E+00
.2230038023239400E-01
-.1817785168087260E+00
.9830866003909881E+00
-.9830866003966130E+00
-.2230038011332300E-01
-.1817785167929150E+00
-.1817785167104710E+00
.9830866004096310E+00
.2230038021145500E-01
-.2230038027285400E-01
-.1817785168366860E+00
-.9830866003849000E+00
-.1817785169360100E+00
-.9830866003681790E+00
-.2230038020037600E-01
-.2230038017080000E-01
.1817785168418750E+00
.9830866003862559E+00
.1817785167109790E+00
-.9830866004090441E+00
.2230038023321200E-01
.2230038021255800E-01
.1817785168040810E+00
-.9830866003922970E+00
.1817785169343840E+00
.9830866003675029E+00
-.2230038024343100E-01
.9830866003916290E+00
-.2230038033237200E-01
.1817785167929960E+00
-.9830866003805701E+00
.2230038033786500E-01
.1817785168521280E+00
-.8979519869718749E+00
.3766956030353650E+00
.2275580184196640E+00
.3766956029275280E+00
.2275580183392060E+00
-.8979519870375030E+00
.8979519869860531E+00
-.3766956030285690E+00
.2275580183749660E+00
.2275580183055540E+00
-.8979519870419040E+00
.3766956029373660E+00
-.3766956028752610E+00
.2275580184552540E+00
.8979519870300200E+00
.2275580184865670E+00
.8979519869904800E+00
-.3766956029506000E+00
-.3766956029825110E+00
-.2275580183689100E+00
-.8979519870069100E+00
-.2275580182809390E+00
.8979519870547670E+00
.3766956029215730E+00
.3766956029314370E+00
-.2275580184255800E+00
.8979519870139741E+00
-.2275580185113490E+00
-.8979519870023480E+00
-.3766956029073390E+00
-.8979519870721940E+00
-.3766956028306370E+00
-.2275580183627070E+00
.8979519870578190E+00
.3766956028230510E+00
-.2275580184319890E+00
-.1713301512452210E+00
.4597861949530550E+00
-.8713453013615680E+00
.4597861948431170E+00
-.8713453014146491E+00
-.1713301512702920E+00
.1713301511912190E+00
-.4597861949823340E+00
-.8713453013567360E+00
-.8713453013647540E+00
-.1713301511629810E+00
.4597861949776620E+00
-.4597861950424320E+00
-.8713453013037380E+00
.1713301512994720E+00
-.8713453013534070E+00
.1713301513627270E+00
-.4597861949247340E+00
-.4597861948552020E+00
.8713453014084100E+00
-.1713301512695920E+00
.8713453013928349E+00
.1713301511781830E+00
.4597861949187800E+00
.4597861950544120E+00
.8713453013090380E+00
.1713301512403680E+00
.8713453013254860E+00
-.1713301513773550E+00
-.4597861949721960E+00
-.1713301512966100E+00
-.4597861949130030E+00
.8713453013725970E+00
.1713301513507360E+00
.4597861949429830E+00
.8713453013461350E+00
-.3971917022972230E+00
-.5480955906492260E+00
-.7360910100912190E+00
-.5480955907789020E+00
-.7360910100565570E+00
-.3971917021825150E+00
.3971917022502210E+00
.5480955906252050E+00
-.7360910101344670E+00
-.7360910101747640E+00
-.3971917021370830E+00
-.5480955906530750E+00
.5480955906102120E+00
-.7360910101691310E+00
.3971917022066690E+00
-.7360910100491940E+00
.3971917023058890E+00
.5480955906993850E+00
.5480955907525290E+00
.7360910100441170E+00
-.3971917022419620E+00
.7360910101399250E+00
.3971917021196020E+00
-.5480955907125310E+00
-.5480955905843859E+00
.7360910101826250E+00
.3971917022173000E+00
.7360910100837820E+00
-.3971917022879800E+00
.5480955906659120E+00
-.3971917021252600E+00
.5480955907404190E+00
.7360910101161060E+00
.3971917021713860E+00
-.5480955907162950E+00
.7360910101091790E+00
.3794747255349560E+00
.6962772780944900E+00
.6092592918368150E+00
.6962772782104411E+00
.6092592917871140E+00
.3794747254020010E+00
-.3794747254955760E+00
-.6962772780741610E+00
.6092592918845759E+00
.6092592919259530E+00
.3794747253762130E+00
.6962772781030080E+00
-.6962772780710560E+00
.6092592919338879E+00
-.3794747254221020E+00
.6092592917959100E+00
-.3794747255155420E+00
-.6962772781408640E+00
-.6962772781859060E+00
-.6092592917748491E+00
.3794747254667130E+00
-.6092592918828780E+00
-.3794747253530890E+00
.6962772781533030E+00
.6962772780465480E+00
-.6092592919465890E+00
-.3794747254466760E+00
-.6092592918387369E+00
.3794747254930950E+00
-.6962772781156230E+00
.3794747253362900E+00
-.6962772781815950E+00
-.6092592918610080E+00
-.3794747253752370E+00
.6962772781612160E+00
-.6092592918600390E+00
-.6787014464703280E+00
.7297642134790810E+00
.8251387328409700E-01
.7297642133897720E+00
.8251387317923400E-01
-.6787014465791040E+00
.6787014464747720E+00
-.7297642134772200E+00
.8251387326399499E-01
.8251387321767099E-01
-.6787014465525470E+00
.7297642134101250E+00
-.7297642133709740E+00
.8251387336840200E-01
.6787014465763180E+00
.8251387332689200E-01
.6787014465346920E+00
-.7297642134143810E+00
-.7297642134312839E+00
-.8251387320173600E-01
-.6787014465317330E+00
-.8251387317169399E-01
.6787014465773989E+00
.7297642133922100E+00
.7297642134127970E+00
-.8251387334668001E-01
.6787014465339900E+00
-.8251387337365500E-01
-.6787014465583360E+00
-.7297642133871040E+00
-.6787014466415410E+00
-.7297642133248270E+00
-.8251387324006100E-01
.6787014466370160E+00
.7297642133213440E+00
-.8251387330807500E-01
]';

x=xyz(3*(j-1)+1);
y=xyz(3*(j-1)+2);
z=xyz(3*(j-1)+3);

th=acos(z);
ph=atan2(y,x);

a=ones(size(th))*(4*pi/84);




function [a,th,ph] = design3_144_16;
% 144 samples, N=8.

j=1:144;

xyz=[
.9383118258138560E+00
-.1750792557749200E+00
-.2981915017822760E+00
-.1751096322456290E+00
-.2982825311210240E+00
.9382772235980340E+00
-.9383116523013461E+00
.1751477614500080E+00
-.2981518150449020E+00
-.2981827578157150E+00
.9383270575537280E+00
-.1750125024219040E+00
.1750977124101310E+00
-.2980583478457380E+00
-.9383506873169580E+00
-.2981854777577620E+00
-.9383236127415390E+00
.1750263369497320E+00
.1751212256614090E+00
.2980709997422250E+00
.9383422805328110E+00
.2981590222823750E+00
-.9382974848874340E+00
-.1752113788700180E+00
-.1751366381351110E+00
.2982885002265250E+00
-.9382702854803310E+00
.2981750565054620E+00
.9382926280748330E+00
.1752101018160420E+00
.9383097216767580E+00
.1750911370548140E+00
.2981911466354040E+00
-.9383070207140820E+00
-.1751442959881740E+00
.2981684263322820E+00
.3183193898656830E+00
-.1895522954118680E+00
.9288394335619220E+00
-.1894661062614570E+00
.9288339463361680E+00
.3183867062421130E+00
-.3182933144730710E+00
.1893628596173800E+00
.9288870078536330E+00
.9288529435535660E+00
.3183507003489590E+00
-.1894334733863170E+00
.1894416073975330E+00
.9288927988957520E+00
-.3182295485121640E+00
.9288662644063450E+00
-.3183138373071290E+00
.1894301027466670E+00
.1894518259149400E+00
-.9288871565521020E+00
.3182399347191460E+00
-.9288657503320540E+00
-.3182891226867960E+00
-.1894741466251780E+00
-.1894810419822530E+00
-.9288341329001750E+00
-.3183772735119440E+00
-.9288638749080860E+00
.3182773954415380E+00
.1895030380803610E+00
.3182754841245910E+00
.1895728238082200E+00
-.9288502897015400E+00
-.3183459025831120E+00
-.1893534180173150E+00
-.9288709110493790E+00
.4152709071162880E+00
.6265468605244530E+00
.6595370385882560E+00
.6266126549472570E+00
.6594514158910070E+00
.4153076097777360E+00
-.4152418281129630E+00
-.6266763943801670E+00
.6594322716641020E+00
.6594942179223080E+00
.4151968471621200E+00
.6266410093775210E+00
-.6266189964270690E+00
.6595218123324770E+00
-.4151862381804330E+00
.6594787856877939E+00
-.4151922150229020E+00
-.6266603193215040E+00
-.6266022331854350E+00
-.6595287758101400E+00
.4152004759696260E+00
-.6594726936833410E+00
-.4153261780732930E+00
.6265779537240910E+00
.6266060528732360E+00
-.6594438365947900E+00
-.4153296051087230E+00
-.6594986338231030E+00
.4153157815166040E+00
-.6265575421369630E+00
.4152509631584860E+00
-.6265428543902710E+00
-.6595534013318720E+00
-.4152672330732850E+00
.6266741585574390E+00
-.6594183983875370E+00
.8147686975402800E-01
.8847674930322230E+00
.4588551001880220E+00
.8848021501705901E+00
.4587806295976860E+00
.8151986849505800E-01
-.8148097265168000E-01
-.8848439651039500E+00
.4587068873636580E+00
.4587780511560210E+00
.8139667888042000E-01
.8848148283368230E+00
-.8848095158928860E+00
.4587845170278200E+00
-.8141798032957800E-01
.4587323275728680E+00
-.8138617295209800E-01
-.8848395009784490E+00
-.8848064690255750E+00
-.4587847998886890E+00
.8144949208920500E-01
-.4587707681467430E+00
-.8156762447812400E-01
.8848028621851550E+00
.8848211768135870E+00
-.4587411019232240E+00
-.8153579869288199E-01
-.4588108991167440E+00
.8157388735636099E-01
-.8847814757064350E+00
.8147060004176100E-01
-.8847779034947540E+00
-.4588361393964780E+00
-.8149754501781800E-01
.8848501095969900E+00
-.4586920902983440E+00
-.7225816121467720E+00
.6911694469079300E+00
-.1267317830534700E-01
.6911462318877840E+00
-.1272247709073500E-01
-.7226029509516230E+00
.7225897391740940E+00
-.6911572322235680E+00
-.1287436155202900E-01
-.1271999109003300E-01
-.7226498291394290E+00
.6910972625263569E+00
-.6911640879369000E+00
-.1283280970189800E-01
.7225839207604250E+00
-.1274089462228200E-01
.7226581266795230E+00
-.6910882009904870E+00
-.6911848254516650E+00
.1280693240541800E-01
-.7225645435168510E+00
.1278690708865000E-01
.7225094351193579E+00
.6912428018942500E+00
.6911667589030220E+00
.1267926954320300E-01
.7225840764307940E+00
.1279840273451600E-01
-.7225177746588930E+00
-.6912338722815931E+00
-.7225871989731150E+00
-.6911634956048890E+00
.1267920436853000E-01
.7225783528006580E+00
.6911703359443890E+00
.1280979212978900E-01
.5601175739954590E+00
.8068680228904130E+00
.1877026822886580E+00
.8068834787163790E+00
.1875714426539700E+00
.5601392734626480E+00
-.5601340935408990E+00
-.8068916312063850E+00
.1875518401461700E+00
.1876521312373620E+00
.5600251494167630E+00
.8069439321680340E+00
-.8068854415129990E+00
.1876857418815800E+00
-.5600981579765580E+00
.1876302229010670E+00
-.5600047208391949E+00
-.8069632036790220E+00
-.8068746775941580E+00
-.1876975169586680E+00
.5601097185238560E+00
-.1876148088020380E+00
-.5602157603217920E+00
.8068202931293010E+00
.8068923207022480E+00
-.1875761333114300E+00
-.5601249655243670E+00
-.1876364876816170E+00
.5602249742367100E+00
-.8068088539003420E+00
.5600918210887200E+00
-.8068808651992270E+00
-.1877243226778800E+00
-.5601293840974760E+00
.8068962450831860E+00
-.1875460549871360E+00
-.9948563422103200E-01
-.3588951295179950E+00
-.9280608248341810E+00
-.3590507942888110E+00
-.9279946080877720E+00
-.9954162185034500E-01
.9943438966061501E-01
.3591437619459990E+00
-.9279701290494740E+00
-.9280190267200989E+00
-.9942019096838001E-01
-.3590213248169130E+00
.3589908155319930E+00
-.9280357484444770E+00
.9937426212442400E-01
-.9280072072034909E+00
.9942025966856401E-01
.3590518560679110E+00
.3590029825622480E+00
.9280313484672880E+00
-.9937139816565700E-01
.9280179389220590E+00
.9951094937970199E-01
-.3589989916314580E+00
-.3590428637423850E+00
.9279961920762100E+00
.9955545935668900E-01
.9280136656320840E+00
-.9948954910509600E-01
.3590159690305810E+00
-.9945187531254500E-01
.3589267513480540E+00
.9280522138670591E+00
.9946550331739699E-01
-.3591200632919870E+00
.9279759661709870E+00
.7878331994376070E+00
.5574500823251660E+00
-.2618554096816970E+00
.5574053886878521E+00
-.2619772920486170E+00
.7878243021845780E+00
-.7878614778767180E+00
-.5573641116878389E+00
-.2619533139327300E+00
-.2618610280706080E+00
.7878026576023160E+00
.5574906049903740E+00
-.5574272044780030E+00
-.2618353048552930E+00
-.7878560686059189E+00
-.2618500918686550E+00
-.7878041469245110E+00
-.5574936371627220E+00
-.5573980474810630E+00
.2618066241900950E+00
.7878862279507650E+00
.2619281468060600E+00
-.7878933745001880E+00
.5573308499710470E+00
.5573998343635920E+00
.2619357785378840E+00
-.7878420352920970E+00
.2619095353283640E+00
.7879083713373940E+00
-.5573183952468600E+00
.7878589677335660E+00
-.5574443214494930E+00
.2617901362647470E+00
-.7878560239272930E+00
.5573693294883240E+00
.2619586152567080E+00
-.5072827321686140E+00
-.7170499460473529E+00
-.4780205063771150E+00
-.7170643140017600E+00
-.4779062710060660E+00
-.5073700481091310E+00
.5073317531927670E+00
.7171162628030801E+00
-.4778689755839950E+00
-.4778916169164080E+00
-.5072575026701600E+00
-.7171536993321960E+00
.7171087443614590E+00
-.4779885176598600E+00
.5072297563685140E+00
-.4779136769269750E+00
.5072353404128420E+00
.7171546742805260E+00
.7171036377589221E+00
.4779420929439370E+00
-.5072807196270021E+00
.4779497913306490E+00
.5073623117813870E+00
-.7170407809538000E+00
-.7170736052366210E+00
.4778899243873540E+00
.5073723138307850E+00
.4779668855044820E+00
-.5073968954290570E+00
.7170049141185160E+00
-.5072894944901550E+00
.7170398743210130E+00
.4780284378712520E+00
.5073429733358930E+00
-.7171476166924809E+00
.4778100075123900E+00
-.4697053900856580E+00
-.3362487640635100E+00
.8162803533040850E+00
-.3361804588591880E+00
.8163540175197370E+00
-.4696262526314000E+00
.4697292672795090E+00
.3360874275716510E+00
.8163330548797630E+00
.8162993201022140E+00
-.4696884809882010E+00
-.3362263386881820E+00
.3361661885920780E+00
.8162610446461910E+00
.4697980423975660E+00
.8163081878411480E+00
.4696844879905110E+00
.3362103868184210E+00
.3361611964247630E+00
-.8162545204851160E+00
-.4698129498065010E+00
-.8163147475476900E+00
.4697491969062010E+00
-.3361040388661380E+00
-.3361667115393140E+00
-.8163550823770680E+00
.4696342422885870E+00
-.8163020291364350E+00
-.4697523383167870E+00
.3361305369549900E+00
-.4697259147648690E+00
.3362543092749910E+00
-.8162662583326020E+00
.4697157090205860E+00
-.3360825711370180E+00
-.8163428557151829E+00
.2209757831175440E+00
.5619818996413199E+00
-.7970859726222270E+00
.5618985433861100E+00
-.7971884426164270E+00
.2208180560990520E+00
-.2209098087123600E+00
-.5618199356383180E+00
-.7972184244822900E+00
-.7971433029261999E+00
.2209065606243460E+00
.5619277943588750E+00
-.5619110464580350E+00
-.7971135607042630E+00
-.2210564344456110E+00
-.7971666081668140E+00
-.2209214541641100E+00
-.5618888768376120E+00
-.5619031892145560E+00
.7971171958991410E+00
.2210632985196790E+00
.7971490712061960E+00
-.2210199177081820E+00
.5618750343727400E+00
.5618715451873800E+00
.7971902229922720E+00
-.2208803184402730E+00
.7971513117794931E+00
.2209665013294830E+00
-.5618928647157230E+00
.2209896742277390E+00
-.5619592284389200E+00
.7970981052900900E+00
-.2209345147362070E+00
.5618214796441769E+00
.7972104899013210E+00
-.2558632109166300E-01
.9914006599926770E+00
-.1283357765359230E+00
.9913910231541920E+00
-.1284105096544480E+00
-.2558476538037500E-01
.2555531861480000E-01
-.9913788670530650E+00
-.1285101850091180E+00
-.1284273557345780E+00
-.2568716764003100E-01
.9913861930235141E+00
-.9913884123582900E+00
-.1284322897283740E+00
.2557660643704000E-01
-.1284710461501210E+00
.2569665752758400E-01
-.9913802863144920E+00
-.9913887704923130E+00
.1284336117570290E+00
-.2555607780520200E-01
.1284346430688090E+00
.2546732329508000E-01
.9913909208299070E+00
.9913860541495390E+00
.1284483459343360E+00
.2558738096289900E-01
.1283929891583590E+00
-.2544975216483000E-01
-.9913967674194480E+00
-.2558970505166500E-01
-.9913982207318930E+00
.1283539439402070E+00
.2557174693595500E-01
.9913764768665120E+00
.1285253549864190E+00
]';

x=xyz(3*(j-1)+1);
y=xyz(3*(j-1)+2);
z=xyz(3*(j-1)+3);

th=acos(z);
ph=atan2(y,x);

a=ones(size(th))*(4*pi/144);



function [a,th,ph] = design3_240_21;
% 240 samples, N=10.

j=1:240;

xyz=[
.8926535357627230E+00
.4125340536573610E+00
-.1816186104542530E+00
.4125340534250320E+00
-.1816186106417820E+00
.8926535358319380E+00
-.8926535358064071E+00
-.4125340536278530E+00
-.1816186103065750E+00
-.1816186106138490E+00
.8926535357404750E+00
.4125340536352400E+00
-.4125340534774350E+00
-.1816186104226540E+00
-.8926535358523040E+00
-.1816186104513840E+00
-.8926535357628120E+00
-.4125340536584320E+00
-.4125340533170900E+00
.1816186106118270E+00
.8926535358879180E+00
.1816186104001360E+00
-.8926535358123000E+00
.4125340535739110E+00
.4125340533279960E+00
.1816186104204000E+00
-.8926535359218249E+00
.1816186105807890E+00
.8926535358109040E+00
-.4125340534973990E+00
.8926535358676440E+00
-.4125340534725580E+00
.1816186103583390E+00
-.8926535358550640E+00
.4125340535351600E+00
.1816186102779710E+00
-.2920937425934330E+00
-.2957670279931700E+00
.9095070701703470E+00
-.2957670280268870E+00
.9095070700892600E+00
-.2920937428117760E+00
.2920937424478640E+00
.2957670280397130E+00
.9095070702019620E+00
.9095070701476120E+00
-.2920937429267210E+00
-.2957670277339340E+00
.2957670281453960E+00
.9095070700844410E+00
.2920937427067830E+00
.9095070701888540E+00
.2920937426892070E+00
.2957670278416750E+00
.2957670279073110E+00
-.9095070701484190E+00
-.2920937427486510E+00
-.9095070701012210E+00
.2920937431592720E+00
-.2957670276469270E+00
-.2957670278353330E+00
-.9095070700472930E+00
.2920937431364140E+00
-.9095070702185910E+00
-.2920937427217760E+00
.2957670277180690E+00
-.2920937425408960E+00
.2957670277931470E+00
-.9095070702522660E+00
.2920937428619380E+00
-.2957670277476140E+00
-.9095070701639690E+00
-.5752257180381920E+00
.2412057282507800E-01
.8176390225974030E+00
.2412057278614400E-01
.8176390225112380E+00
-.5752257181623009E+00
.5752257181164780E+00
-.2412057297921300E-01
.8176390225377810E+00
.8176390225560030E+00
-.5752257181034800E+00
.2412057267146900E-01
-.2412057304150300E-01
.8176390224407570E+00
.5752257182517770E+00
.8176390224583791E+00
.5752257182291181E+00
-.2412057298452600E-01
-.2412057281823900E-01
-.8176390225812600E+00
-.5752257180614240E+00
-.8176390225435780E+00
.5752257181238820E+00
.2412057260611100E-01
.2412057271295000E-01
-.8176390225272960E+00
.5752257181425460E+00
-.8176390226004950E+00
-.5752257180351740E+00
-.2412057279222800E-01
-.5752257179254689E+00
-.2412057271105200E-01
-.8176390226800700E+00
.5752257179082300E+00
.2412057259415500E-01
-.8176390226956460E+00
-.1288331617248000E+00
.5224764072024000E-01
.9902889479738530E+00
.5224764069440900E-01
.9902889479588950E+00
-.1288331618502510E+00
.1288331618403250E+00
-.5224764032003800E-01
.9902889479799381E+00
.9902889479497170E+00
-.1288331619247960E+00
.5224764068455800E-01
-.5224764038851000E-01
.9902889479675810E+00
.1288331619075380E+00
.9902889479777300E+00
.1288331618780010E+00
-.5224764026899200E-01
-.5224764039040900E-01
-.9902889479621900E+00
-.1288331619482090E+00
-.9902889479606261E+00
.1288331618966490E+00
.5224764054718700E-01
.5224764052780800E-01
-.9902889479532510E+00
.1288331619612000E+00
-.9902889479708680E+00
-.1288331619362050E+00
-.5224764025552600E-01
-.1288331617904780E+00
-.5224764033764300E-01
-.9902889479854941E+00
.1288331618574160E+00
.5224764055154500E-01
-.9902889479655000E+00
.7180063860347500E+00
.6574468762559930E+00
-.2285397875962860E+00
.6574468762867370E+00
-.2285397878319220E+00
.7180063859315960E+00
-.7180063861094420E+00
-.6574468761714340E+00
-.2285397876048770E+00
-.2285397877372190E+00
.7180063859474221E+00
.6574468763023740E+00
-.6574468762410210E+00
-.2285397877138000E+00
-.7180063860110540E+00
-.2285397876789970E+00
-.7180063860313590E+00
-.6574468762309450E+00
-.6574468763611850E+00
.2285397878605490E+00
.7180063858543150E+00
.2285397877030650E+00
-.7180063858573850E+00
.6574468764125770E+00
.6574468763044540E+00
.2285397878740170E+00
-.7180063859019750E+00
.2285397877849670E+00
.7180063858138530E+00
-.6574468764316480E+00
.7180063858807600E+00
-.6574468763634850E+00
.2285397877708510E+00
-.7180063858910180E+00
.6574468763715580E+00
.2285397877154010E+00
.8631764731178030E+00
.4681818166531380E+00
.1890295289400010E+00
.4681818164384860E+00
.1890295291974920E+00
.8631764731778400E+00
-.8631764731944460E+00
-.4681818165764200E+00
.1890295287800330E+00
.1890295291255270E+00
.8631764730643890E+00
.4681818166767080E+00
-.4681818163926710E+00
.1890295288974430E+00
-.8631764732683980E+00
.1890295287921740E+00
-.8631764731436880E+00
-.4681818166651000E+00
-.4681818164112130E+00
-.1890295291281380E+00
.8631764732078210E+00
-.1890295288978520E+00
-.8631764730897200E+00
.4681818167219310E+00
.4681818165088670E+00
-.1890295289305550E+00
-.8631764731981230E+00
-.1890295290018230E+00
.8631764731066590E+00
-.4681818166487220E+00
.8631764731352290E+00
-.4681818166486420E+00
-.1890295288715610E+00
-.8631764731233340E+00
.4681818166987620E+00
-.1890295288017440E+00
.7726328568471330E+00
-.5170594506955900E+00
.3683585114621520E+00
-.5170594505671320E+00
.3683585115855150E+00
.7726328568742860E+00
-.7726328568060810E+00
.5170594506473910E+00
.3683585116159150E+00
.3683585116480010E+00
.7726328568060540E+00
-.5170594506245730E+00
.5170594504940070E+00
.3683585118165880E+00
-.7726328568130560E+00
.3683585117204960E+00
-.7726328568024760E+00
.5170594505782730E+00
.5170594505834450E+00
-.3683585114871170E+00
.7726328569102801E+00
-.3683585115673300E+00
-.7726328568594670E+00
-.5170594506022290E+00
-.5170594505023690E+00
-.3683585116659560E+00
-.7726328568792750E+00
-.3683585114698030E+00
.7726328568556510E+00
.5170594506774120E+00
.7726328569347490E+00
.5170594506919191E+00
-.3683585112835310E+00
-.7726328569274850E+00
-.5170594506337780E+00
-.3683585113803780E+00
-.8478192319146480E+00
-.6632577590016700E-01
-.5261211281130020E+00
-.6632577591363099E-01
-.5261211282576860E+00
-.8478192318238090E+00
.8478192318830180E+00
.6632577581985200E-01
-.5261211281740970E+00
-.5261211283487620E+00
-.8478192317669569E+00
-.6632577591791000E-01
.6632577584612000E-01
-.5261211284070980E+00
.8478192317363720E+00
-.5261211284592400E+00
.8478192317090800E+00
.6632577578136600E-01
.6632577594578500E-01
.5261211283443800E+00
-.8478192317674960E+00
.5261211284495320E+00
.8478192317006920E+00
-.6632577596561300E-01
-.6632577587721100E-01
.5261211283063880E+00
.8478192317964360E+00
.5261211285046690E+00
-.8478192316652130E+00
.6632577598176000E-01
-.8478192318217250E+00
.6632577594100500E-01
.5261211282575941E+00
.8478192318502640E+00
-.6632577599665500E-01
.5261211282045900E+00
.9805743229230000E-02
.9429838158425931E+00
.3326941094438920E+00
.9429838158089230E+00
.3326941095397480E+00
.9805743214949999E-02
-.9805743379690000E-02
-.9429838157872910E+00
.3326941095962070E+00
.3326941092265540E+00
.9805743204272001E-02
.9429838159195320E+00
-.9429838157740400E+00
.3326941096356470E+00
-.9805743315804000E-02
.3326941093979960E+00
-.9805743298910000E-02
-.9429838158580620E+00
-.9429838157761140E+00
-.3326941096300980E+00
.9805743304667000E-02
-.3326941093190270E+00
-.9805743188507000E-02
.9429838158870700E+00
.9429838157750819E+00
-.3326941096351990E+00
-.9805743230762999E-02
-.3326941094557650E+00
.9805743389761999E-02
-.9429838158367350E+00
.9805743301140001E-02
-.9429838157525240E+00
-.3326941096970650E+00
-.9805743287713000E-02
.9429838157913790E+00
-.3326941095873310E+00
.7855992483711520E+00
-.4051569453122690E+00
-.4676341204658960E+00
-.4051569449321250E+00
-.4676341206498590E+00
.7855992484576980E+00
-.7855992482017899E+00
.4051569454340510E+00
-.4676341206449040E+00
-.4676341206112420E+00
.7855992483346230E+00
-.4051569452153390E+00
.4051569451364230E+00
-.4676341208682010E+00
-.7855992482223660E+00
-.4676341208118040E+00
-.7855992481456090E+00
.4051569453503470E+00
.4051569448419850E+00
.4676341208613320E+00
.7855992483783050E+00
.4676341207867260E+00
-.7855992482498571E+00
-.4051569451771560E+00
-.4051569449996430E+00
.4676341208710980E+00
-.7855992482911820E+00
.4676341208937130E+00
.7855992482342400E+00
.4051569450839530E+00
.7855992483133410E+00
.4051569451171040E+00
.4676341207321060E+00
-.7855992482810999E+00
-.4051569451973700E+00
.4676341207167270E+00
-.7373319991314921E+00
.6208515010137640E+00
-.2662422519918900E+00
.6208515009491860E+00
-.2662422521548950E+00
-.7373319991270100E+00
.7373319990600610E+00
-.6208515010887370E+00
-.2662422520148830E+00
-.2662422519486310E+00
-.7373319991032550E+00
.6208515010658500E+00
-.6208515010792210E+00
-.2662422522338000E+00
.7373319989890250E+00
-.2662422520116240E+00
.7373319989962220E+00
-.6208515011659510E+00
-.6208515010721239E+00
.2662422522225600E+00
-.7373319989990600E+00
.2662422521138640E+00
.7373319988329740E+00
.6208515013159830E+00
.6208515011873870E+00
.2662422523283740E+00
.7373319988637970E+00
.2662422519322500E+00
-.7373319989389900E+00
-.6208515012679590E+00
-.7373319989479430E+00
-.6208515011832970E+00
.2662422521048790E+00
.7373319988350070E+00
.6208515013057860E+00
.2662422521320100E+00
.7268714691656590E+00
-.2748828235042800E-01
-.6862231864680610E+00
-.2748828218275500E-01
-.6862231864483250E+00
.7268714691906329E+00
-.7268714691729310E+00
.2748828237188500E-01
-.6862231864594990E+00
-.6862231864497120E+00
.7268714691854060E+00
-.2748828228634100E-01
.2748828235160700E-01
-.6862231864911200E+00
-.7268714691438450E+00
-.6862231865456220E+00
-.7268714690897941E+00
.2748828242028100E-01
.2748828226683600E-01
.6862231864703350E+00
.7268714691666740E+00
.6862231866611830E+00
-.7268714689834220E+00
-.2748828234818500E-01
-.2748828225102900E-01
.6862231865230920E+00
-.7268714691174650E+00
.6862231866091120E+00
.7268714690334980E+00
.2748828232394800E-01
.7268714690701070E+00
.2748828233555000E-01
.6862231865698690E+00
-.7268714690801830E+00
-.2748828230971600E-01
.6862231865602320E+00
.6653633857205150E+00
.5808602677392710E+00
.4689274083527160E+00
.5808602675770870E+00
.4689274084886380E+00
.6653633857663080E+00
-.6653633856773800E+00
-.5808602677195750E+00
.4689274084383180E+00
.4689274083407830E+00
.6653633858218631E+00
.5808602676328130E+00
-.5808602675284530E+00
.4689274086788320E+00
-.6653633856747230E+00
.4689274083726140E+00
-.6653633856988030E+00
-.5808602677480780E+00
-.5808602676408769E+00
-.4689274085527620E+00
.6653633856654270E+00
-.4689274084683360E+00
-.6653633858479470E+00
.5808602674999610E+00
.5808602673867520E+00
-.4689274086545190E+00
-.6653633858155630E+00
-.4689274083756990E+00
.6653633856513560E+00
-.5808602677999380E+00
.6653633856518190E+00
-.5808602677912120E+00
-.4689274083858500E+00
-.6653633857517340E+00
.5808602675480170E+00
-.4689274085453260E+00
-.5801253673053040E+00
-.7790995979244340E+00
.2376097109187070E+00
-.7790995980535180E+00
.2376097109099340E+00
-.5801253671355390E+00
.5801253671868080E+00
.7790995979777320E+00
.2376097110332580E+00
.2376097106959320E+00
-.5801253672761100E+00
-.7790995980141140E+00
.7790995980647319E+00
.2376097111473200E+00
.5801253670232500E+00
.2376097108192850E+00
.5801253670474260E+00
.7790995981467740E+00
.7790995981702240E+00
-.2376097108496420E+00
-.5801253670034990E+00
-.2376097108118020E+00
.5801253671572560E+00
-.7790995980672760E+00
-.7790995980749610E+00
-.2376097110451280E+00
.5801253670513690E+00
-.2376097106092530E+00
-.5801253670223590E+00
.7790995982294950E+00
-.5801253670900940E+00
.7790995981519661E+00
-.2376097106980860E+00
.5801253672184110E+00
-.7790995979667160E+00
-.2376097109922150E+00
.9586680253602000E+00
.1011136059005390E+00
-.2659542363899560E+00
.1011136058898930E+00
-.2659542364771990E+00
.9586680253371200E+00
-.9586680253264100E+00
-.1011136060954320E+00
-.2659542364376600E+00
-.2659542366341790E+00
.9586680252945550E+00
.1011136058805580E+00
-.1011136060031710E+00
-.2659542366563170E+00
-.9586680252754820E+00
-.2659542367154550E+00
-.9586680252461620E+00
-.1011136061256020E+00
-.1011136058254380E+00
.2659542364146640E+00
.9586680253612671E+00
.2659542362867390E+00
-.9586680253935830E+00
.1011136058555220E+00
.1011136058024440E+00
.2659542362606640E+00
-.9586680254064150E+00
.2659542365158540E+00
.9586680253225770E+00
-.1011136059261060E+00
.9586680254495000E+00
-.1011136059091010E+00
.2659542360648080E+00
-.9586680254786000E+00
.1011136057894970E+00
.2659542360053860E+00
-.7844318144170850E+00
.2843190250072290E+00
.5512072392025160E+00
.2843190248228480E+00
.5512072393207090E+00
-.7844318144008620E+00
.7844318144434220E+00
-.2843190248881310E+00
.5512072392264670E+00
.5512072394346770E+00
-.7844318142918880E+00
.2843190249025560E+00
-.2843190246401610E+00
.5512072393475040E+00
.7844318144482491E+00
.5512072394083570E+00
.7844318144009980E+00
-.2843190246525460E+00
-.2843190247149400E+00
-.5512072391601370E+00
-.7844318145528040E+00
-.5512072394176490E+00
.7844318144267430E+00
.2843190245635030E+00
.2843190244771060E+00
-.5512072393940670E+00
.7844318144746290E+00
-.5512072392271640E+00
-.7844318145108320E+00
-.2843190247007970E+00
-.7844318146549000E+00
-.2843190247577290E+00
-.5512072389927720E+00
.7844318145421389E+00
.2843190246898840E+00
-.5512072391882400E+00
.1666638785351180E+00
.9794687788666500E+00
.1134198519532850E+00
.9794687788923619E+00
.1134198520112480E+00
.1666638783445640E+00
-.1666638783223350E+00
-.9794687788772219E+00
.1134198521746590E+00
.1134198518526030E+00
.1666638784650920E+00
.9794687788902240E+00
-.9794687789080510E+00
.1134198522332290E+00
-.1666638781012970E+00
.1134198520235320E+00
-.1666638782131650E+00
-.9794687789132980E+00
-.9794687788914180E+00
-.1134198520887550E+00
.1666638782973680E+00
-.1134198519422990E+00
-.1666638783837850E+00
.9794687788936730E+00
.9794687788877920E+00
-.1134198522526510E+00
-.1666638782071420E+00
-.1134198518873330E+00
.1666638784200610E+00
-.9794687788938650E+00
.1666638785133120E+00
-.9794687788588400E+00
-.1134198520527750E+00
-.1666638785259920E+00
.9794687788524030E+00
-.1134198520897270E+00
.9035426353908700E+00
.9900269067959901E-01
.4169042735078650E+00
.9900269051118001E-01
.4169042737536920E+00
.9035426352958970E+00
-.9035426353835330E+00
-.9900269064792300E-01
.4169042735312880E+00
.4169042739582500E+00
.9035426351937680E+00
.9900269058185000E-01
-.9900269041493300E-01
.4169042736997320E+00
-.9035426353313410E+00
.4169042738439640E+00
-.9035426352375170E+00
-.9900269066384500E-01
-.9900269046419199E-01
-.4169042739372540E+00
.9035426352163480E+00
-.4169042742060360E+00
-.9035426351101470E+00
.9900269030157501E-01
.9900269012804400E-01
-.4169042740643800E+00
-.9035426351945230E+00
-.4169042741137440E+00
.9035426351313860E+00
-.9900269049639200E-01
.9035426352792750E+00
-.9900269046710200E-01
-.4169042738001830E+00
-.9035426352343990E+00
.9900269024582900E-01
-.4169042739499880E+00
.2787624045360920E+00
.3493121855370630E+00
-.8945795206981750E+00
.3493121855860560E+00
-.8945795206085150E+00
.2787624047624310E+00
-.2787624045405250E+00
-.3493121855034730E+00
-.8945795207099100E+00
-.8945795207341440E+00
.2787624047279170E+00
.3493121852918660E+00
-.3493121854667010E+00
-.8945795206777230E+00
-.2787624046898960E+00
-.8945795207888640E+00
-.2787624046586770E+00
-.3493121852069840E+00
-.3493121855510410E+00
.8945795206827980E+00
.2787624045679230E+00
.8945795207852190E+00
-.2787624046804690E+00
.3493121851989290E+00
.3493121855496230E+00
.8945795206792300E+00
-.2787624045811490E+00
.8945795207818050E+00
.2787624045559080E+00
-.3493121853070750E+00
.2787624044379500E+00
-.3493121855065000E+00
.8945795207406920E+00
-.2787624044432590E+00
.3493121854287870E+00
.8945795207693820E+00
.5558962301794150E+00
-.6768332117366710E+00
.4825724658147600E+00
-.6768332116815670E+00
.4825724660401160E+00
.5558962300508760E+00
-.5558962303148920E+00
.6768332115229870E+00
.4825724659584010E+00
.4825724659102830E+00
.5558962301646720E+00
-.6768332116806730E+00
.6768332114576920E+00
.4825724660928950E+00
-.5558962302776390E+00
.4825724659029810E+00
-.5558962303679090E+00
.6768332115189570E+00
.6768332116355920E+00
-.4825724660719810E+00
.5558962300791910E+00
-.4825724661505860E+00
-.5558962302300841E+00
-.6768332114556160E+00
-.6768332114382860E+00
-.4825724663277370E+00
-.5558962300974000E+00
-.4825724659723730E+00
.5558962302677700E+00
.6768332115517270E+00
.5558962301926910E+00
.6768332115894530E+00
-.4825724660059490E+00
-.5558962301943380E+00
-.6768332114555370E+00
-.4825724661918750E+00
]';

x=xyz(3*(j-1)+1);
y=xyz(3*(j-1)+2);
z=xyz(3*(j-1)+3);

th=acos(z);
ph=atan2(y,x);

a=ones(size(th))*(4*pi/240);


