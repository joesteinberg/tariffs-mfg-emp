**DST project

*********************Aggregating 6 digit HS codes to 2 digit ISIC industry codes

{gen twodigit=1 if commoditycode>=10000 & commoditycode<=19999

replace twodigit=2 if commoditycode>=20000 & commoditycode<=29999

replace twodigit=3 if commoditycode>=30000 & commoditycode<=39999

replace twodigit=4 if commoditycode>=40000 & commoditycode<=49999

replace twodigit=5 if commoditycode>=50000 & commoditycode<=59999

replace twodigit=6 if commoditycode>=60000 & commoditycode<=69999

replace twodigit=7 if commoditycode>=70000 & commoditycode<=79999

replace twodigit=8 if commoditycode>=80000 & commoditycode<=89999

replace twodigit=9 if commoditycode>=90000 & commoditycode<=99999

**

replace twodigit=10 if commoditycode>=100000 & commoditycode<=109999

replace twodigit=11 if commoditycode>=110000 & commoditycode<=119999

replace twodigit=12 if commoditycode>=120000 & commoditycode<=129999

replace twodigit=13 if commoditycode>=130000 & commoditycode<=139999

replace twodigit=14 if commoditycode>=140000 & commoditycode<=149999

replace twodigit=15 if commoditycode>=150000 & commoditycode<=159999

replace twodigit=16 if commoditycode>=160000 & commoditycode<=169999

replace twodigit=17 if commoditycode>=170000 & commoditycode<=179999

replace twodigit=18 if commoditycode>=180000 & commoditycode<=189999

replace twodigit=19 if commoditycode>=190000 & commoditycode<=199999

**

replace twodigit=20 if commoditycode>=200000 & commoditycode<=209999

replace twodigit=21 if commoditycode>=210000 & commoditycode<=219999

replace twodigit=22 if commoditycode>=220000 & commoditycode<=229999

replace twodigit=23 if commoditycode>=230000 & commoditycode<=239999

replace twodigit=24 if commoditycode>=240000 & commoditycode<=249999

replace twodigit=25 if commoditycode>=250000 & commoditycode<=259999

replace twodigit=26 if commoditycode>=260000 & commoditycode<=269999

replace twodigit=27 if commoditycode>=270000 & commoditycode<=279999

replace twodigit=28 if commoditycode>=280000 & commoditycode<=289999

replace twodigit=29 if commoditycode>=290000 & commoditycode<=299999

replace twodigit=30 if commoditycode>=300000 & commoditycode<=309999

**

replace twodigit=31 if commoditycode>=310000 & commoditycode<=319999

replace twodigit=32 if commoditycode>=320000 & commoditycode<=329999

replace twodigit=33 if commoditycode>=330000 & commoditycode<=339999

replace twodigit=34 if commoditycode>=340000 & commoditycode<=349999

replace twodigit=35 if commoditycode>=350000 & commoditycode<=359999

replace twodigit=36 if commoditycode>=360000 & commoditycode<=369999

replace twodigit=37 if commoditycode>=370000 & commoditycode<=379999

replace twodigit=38 if commoditycode>=380000 & commoditycode<=389999

replace twodigit=39 if commoditycode>=390000 & commoditycode<=399999

replace twodigit=40 if commoditycode>=400000 & commoditycode<=409999

**

replace twodigit=41 if commoditycode>=410000 & commoditycode<=419999

replace twodigit=42 if commoditycode>=420000 & commoditycode<=429999

replace twodigit=43 if commoditycode>=430000 & commoditycode<=439999

replace twodigit=44 if commoditycode>=440000 & commoditycode<=449999

replace twodigit=45 if commoditycode>=450000 & commoditycode<=459999

replace twodigit=46 if commoditycode>=460000 & commoditycode<=469999

replace twodigit=47 if commoditycode>=470000 & commoditycode<=479999

replace twodigit=48 if commoditycode>=480000 & commoditycode<=489999

replace twodigit=49 if commoditycode>=490000 & commoditycode<=499999

replace twodigit=50 if commoditycode>=500000 & commoditycode<=509999

**

replace twodigit=51 if commoditycode>=510000 & commoditycode<=519999

replace twodigit=52 if commoditycode>=520000 & commoditycode<=529999

replace twodigit=53 if commoditycode>=530000 & commoditycode<=539999

replace twodigit=54 if commoditycode>=540000 & commoditycode<=549999

replace twodigit=55 if commoditycode>=550000 & commoditycode<=559999

replace twodigit=56 if commoditycode>=560000 & commoditycode<=569999

replace twodigit=57 if commoditycode>=570000 & commoditycode<=579999

replace twodigit=58 if commoditycode>=580000 & commoditycode<=589999

replace twodigit=59 if commoditycode>=590000 & commoditycode<=599999

replace twodigit=60 if commoditycode>=600000 & commoditycode<=609999

**

replace twodigit=61 if commoditycode>=610000 & commoditycode<=619999

replace twodigit=62 if commoditycode>=620000 & commoditycode<=629999

replace twodigit=63 if commoditycode>=630000 & commoditycode<=639999

replace twodigit=64 if commoditycode>=640000 & commoditycode<=649999

replace twodigit=65 if commoditycode>=650000 & commoditycode<=659999

replace twodigit=66 if commoditycode>=660000 & commoditycode<=669999

replace twodigit=67 if commoditycode>=670000 & commoditycode<=679999

replace twodigit=68 if commoditycode>=680000 & commoditycode<=689999

replace twodigit=69 if commoditycode>=690000 & commoditycode<=699999

replace twodigit=70 if commoditycode>=700000 & commoditycode<=709999

**

replace twodigit=71 if commoditycode>=710000 & commoditycode<=719999

replace twodigit=72 if commoditycode>=720000 & commoditycode<=729999

replace twodigit=73 if commoditycode>=730000 & commoditycode<=739999

replace twodigit=74 if commoditycode>=740000 & commoditycode<=749999

replace twodigit=75 if commoditycode>=750000 & commoditycode<=759999

replace twodigit=76 if commoditycode>=760000 & commoditycode<=769999

replace twodigit=77 if commoditycode>=770000 & commoditycode<=779999

replace twodigit=78 if commoditycode>=780000 & commoditycode<=789999

replace twodigit=79 if commoditycode>=790000 & commoditycode<=799999

replace twodigit=80 if commoditycode>=800000 & commoditycode<=809999

**

replace twodigit=71 if commoditycode>=710000 & commoditycode<=719999

replace twodigit=72 if commoditycode>=720000 & commoditycode<=729999

replace twodigit=73 if commoditycode>=730000 & commoditycode<=739999

replace twodigit=74 if commoditycode>=740000 & commoditycode<=749999

replace twodigit=75 if commoditycode>=750000 & commoditycode<=759999

replace twodigit=76 if commoditycode>=760000 & commoditycode<=769999

replace twodigit=77 if commoditycode>=770000 & commoditycode<=779999

replace twodigit=78 if commoditycode>=780000 & commoditycode<=789999

replace twodigit=79 if commoditycode>=790000 & commoditycode<=799999

replace twodigit=80 if commoditycode>=800000 & commoditycode<=809999

**

replace twodigit=81 if commoditycode>=810000 & commoditycode<=819999

replace twodigit=82 if commoditycode>=820000 & commoditycode<=829999

replace twodigit=83 if commoditycode>=830000 & commoditycode<=839999

replace twodigit=84 if commoditycode>=840000 & commoditycode<=849999

replace twodigit=85 if commoditycode>=850000 & commoditycode<=859999

replace twodigit=86 if commoditycode>=860000 & commoditycode<=869999

replace twodigit=87 if commoditycode>=870000 & commoditycode<=879999

replace twodigit=88 if commoditycode>=880000 & commoditycode<=889999

replace twodigit=89 if commoditycode>=890000 & commoditycode<=899999

replace twodigit=90 if commoditycode>=900000 & commoditycode<=909999

**

replace twodigit=91 if commoditycode>=910000 & commoditycode<=919999

replace twodigit=92 if commoditycode>=920000 & commoditycode<=929999

replace twodigit=93 if commoditycode>=930000 & commoditycode<=939999

replace twodigit=94 if commoditycode>=940000 & commoditycode<=949999

replace twodigit=95 if commoditycode>=950000 & commoditycode<=959999

replace twodigit=96 if commoditycode>=960000 & commoditycode<=969999

replace twodigit=97 if commoditycode>=970000 & commoditycode<=979999

replace twodigit=98 if commoditycode>=980000 & commoditycode<=989999

replace twodigit=99 if commoditycode>=990000 & commoditycode<=999999}

*************************************

**Dropping of irrelevant variables

drop aggregatelevel isleafcode tradeflowcode reportercode partnercode qtyunitcode qtyunit qty altqtyunitcode altqtyunit altqty netweightkg grossweightkg ciftradevalueus fobtradevalueus flag

**Alternate

drop period perioddesc aggregatelevel isleafcode tradeflowcode reportercode reporteriso partnercode partneriso ndpartnercode ndpartner ndpartneriso customsproccode customs modeoftransportcode modeoftransport qtyunitcode qtyunit qty altqtyunitcode altqtyunit altqty netweightkg grossweightkg ciftradevalueus fobtradevalueus flag


drop if commoditycode<=20000

***HS 1996 to ISIC Rev 4 (also NIC 2008) ###ONLY MANUFACTURING --> TOTAL 22 sectors except code=33 (repair and installation of electrical equipments)

***************** (1)

gen isic4=10 if twodigit==2 |  twodigit==3 | twodigit==4 | twodigit==5 | twodigit==7 | twodigit==8  | twodigit==9 | twodigit==10 |  ///
twodigit==11 | twodigit==12 | twodigit==14 | twodigit==15 | twodigit==16 | twodigit==17 | twodigit==18 | twodigit==19 | twodigit==20 | twodigit==21  | ///
twodigit==23 //manufacturing of food products 

replace isic4=10 if  commoditycode>=410110 & commoditycode<=410310 

replace isic4=10 if commoditycode==350110 | commoditycode==350190 | commoditycode==350211 | commoditycode==350219 |commoditycode==350510  

***************** (2) 

replace isic4=11 if twodigit==22 // manufacturing of beverages

*****************

replace isic4=12 if twodigit==24 // manufacture of tobacco products


***************** (3)

replace isic4=13 if twodigit==50 | twodigit==51 | twodigit==52 | twodigit==53 | twodigit==54 | twodigit==55 | twodigit==56 | twodigit==57 | /// 
twodigit==58 | twodigit==59  | twodigit==60 | twodigit==63 | twodigit==94  // manufacture of textiles

replace isic4=13 if commoditycode==701940 | commoditycode==701951 | commoditycode==701952 | commoditycode==701959

replace isic4=13 if commoditycode==880400 


***************** (4)

replace isic4=14 if twodigit==61 | twodigit==62 | twodigit==65 | twodigit==43 // manufacture of wearing apparel

replace isic4=14 if commoditycode==420310 | commoditycode==420329 | commoditycode==420330 | commoditycode==420340 


***************** (5)

replace isic4=15 if  commoditycode>=410410 & commoditycode<=411510 /// Manufacture of leather and related products

replace isic4=15 if  commoditycode>=420100 & commoditycode<=420299

replace isic4=15 if  commoditycode==420400 | commoditycode==420500 | commoditycode==911390 | commoditycode==960500

replace isic4=15 if twodigit==64 

*****************(6)

replace isic4=16 if twodigit==44 | twodigit==45 | twodigit==46 /// manufacture of wood and of products of wood and cork, except furniture

*****************(7)

replace isic4=17 if twodigit==47 /// manufacture of paper and paper products

replace isic4=17 if commoditycode>=480100 & commoditycode<=482390

replace isic4=17 if commoditycode==590500

*****************(8)

replace isic4=18 if twodigit==49 /// printing and reproduction of recorded media

replace isic4=18 if commoditycode>=482010 & commoditycode<=482090

replace isic4=18 if commoditycode==844250 | commoditycode==852351 | commoditycode==852359 | commoditycode==852380

*****************(9)

replace isic4=19 if commoditycode>=271000 & commoditycode<=271390 /// manufacture of coke and refined petroleum products

replace isic4=19 if commoditycode==270400 | commoditycode==270600

replace isic4=19 if commoditycode==284410 | commoditycode==284420 | commoditycode==284430 | commoditycode==284440

*****************(10)

replace isic4=20 if twodigit==26 | twodigit==29 | twodigit==32 | twodigit==31 | twodigit==33 | twodigit==34 | twodigit==37 /// manufacture of chemicals and chemical products

replace isic4=20 if twodigit==35 & isic4!=15

replace isic4=20 if commoditycode>=270710 & commoditycode<=270820

replace isic4=20 if commoditycode>=280110 & commoditycode<=284390

replace isic4=20 if commoditycode>=284510 & commoditycode<=285300

replace isic4=20 if commoditycode>=380210 & commoditycode<=382590 ///note except 381600 and 382450

replace isic4=20 if commoditycode==440200 | commoditycode==300670 | commoditycode==710410 | commoditycode==710420

replace isic4=20 if commoditycode>=390110 & commoditycode<=391400

replace isic4=20 if commoditycode>=400211 & commoditycode<=400299

replace isic4=20 if commoditycode>=360100 & commoditycode<=360490

replace isic4=20 if commoditycode>=852321 & commoditycode<=852340

replace isic4=20 if commoditycode>=540210 & commoditycode<=540259

replace isic4=20 if commoditycode>=540310 & commoditycode<=540339

replace isic4=20 if commoditycode>=540410 & commoditycode<=540500

replace isic4=20 if commoditycode>=550110 & commoditycode<=550490

replace isic4=20 if commoditycode==151800 | commoditycode==152000

****************(11)

replace isic4=21 if commoditycode==291821 | commoditycode==291822 | commoditycode==291823 | commoditycode==300692 /// manufacture of basic pharmaceutical products and pharmaceutical preparations

replace isic4=21 if commoditycode==292241 | commoditycode==292242 

replace isic4=21 if commoditycode>=292310 & commoditycode<=292410

replace isic4=21 if commoditycode==292422 | commoditycode==292429 | commoditycode==293229 | commoditycode==293311 | commoditycode==293319 | commoditycode==293321 | commoditycode==293430

replace isic4=21 if commoditycode>=293351 & commoditycode<=293369

replace isic4=21 if commoditycode>=293500 & commoditycode<=294190

replace isic4=21 if commoditycode>=300120 & commoditycode<=300660

***************(12)

replace isic4=22 if commoditycode>=391610 & commoditycode<=392690 /// manufacture of rubber and plastics products

replace isic4=22 if commoditycode>=400300 & commoditycode<=401700

replace isic4=22 if commoditycode==590610 | commoditycode==300691 | commoditycode==590691 | commoditycode==590699 | commoditycode==853670

replace isic4=22 if commoditycode==650610 | commoditycode==650691 | commoditycode==854720 | commoditycode==940592 

***************(13)

replace isic4=23 if twodigit==25 | twodigit==68 | twodigit==69 /// manufacture of other non-metallic mineral products

replace isic4=23 if commoditycode==271500 | commoditycode==281810 | commoditycode==380110 | commoditycode==380120 | ///
commoditycode==380130 | commoditycode==380190 | commoditycode==381600 | commoditycode==382450 | commoditycode==854610 | /// 
commoditycode==854620 | commoditycode==854710

replace isic4=23 if twodigit==70 & isic4!=13

***************(14)

replace isic4=24 if twodigit==72 | twodigit==75 | twodigit==78 | twodigit==79 | twodigit==80 | twodigit==81 /// manufacture of basic metals

replace isic4=24 if commoditycode>=730110 & commoditycode<=730799

replace isic4=24 if commoditycode>=740110 & commoditycode<=741220

replace isic4=24 if commoditycode>=760110 & commoditycode<=760900

replace isic4=24 if commoditycode>=710610 & commoditycode<=711100

replace isic4=24 if commoditycode==281820

***************(15)

replace isic4=25 if commoditycode>=730810 & commoditycode<=731450  /// manufacture of fabricated metal products, except machinery and equipment

replace isic4=25 if commoditycode>=731520 & commoditycode<=732090

replace isic4=25 if commoditycode>=732211 & commoditycode<=732290

replace isic4=25 if commoditycode>=732310 & commoditycode<=732690

replace isic4=25 if twodigit==74 & commoditycode!=741700

replace isic4=25 if twodigit==82 | twodigit==83

replace isic4=25 if commoditycode==750810 | commoditycode==750890 | commoditycode==780600 | commoditycode==790700 | commoditycode==800700 | ///
commoditycode==840110 | commoditycode==840410 | commoditycode==840420 | commoditycode==840490 | commoditycode==848710 | commoditycode==940600 

replace isic4=25 if commoditycode>=761010 & commoditycode<=761699

replace isic4=25 if commoditycode>=840140 & commoditycode<=840290

***************(16)

replace isic4=26 if commoditycode==841920 | commoditycode==844312 | commoditycode==852352 | commoditycode==940210 | commoditycode==940290 /// manufacture of computer, electronic and optical products

replace isic4=26 if commoditycode==844331 | commoditycode==844332 | commoditycode==844339 | commoditycode==844399

replace isic4=26 if commoditycode>=846900 & commoditycode<=847350

replace isic4=26 if commoditycode>=851711 & commoditycode<=852290

replace isic4=26 if commoditycode>=852510 & commoditycode<=852990

replace isic4=26 if commoditycode>=853210 & commoditycode<=853400

replace isic4=26 if commoditycode>=854011 & commoditycode<=854290

replace isic4=26 if twodigit==90 & commoditycode!=902300

replace isic4=26 if twodigit==91 & commoditycode!=911390

**************(17)

replace isic4=27 if commoditycode==630110 | commoditycode==732290 | commoditycode==741700 | commoditycode==840310 | /// manufacture of electrical equipments
commoditycode==840390 | commoditycode==841451 | commoditycode==841460 | commoditycode==841911 | ///
commoditycode==841919 | commoditycode==842211 | commoditycode==842219 | commoditycode==845121 | ///
commoditycode==854690 | commoditycode==940599 | commoditycode==900662 

replace isic4=27 if commoditycode>=732111 & commoditycode<=732190

replace isic4=27 if commoditycode>=841810 & commoditycode<=841840

replace isic4=27 if commoditycode>=845011 & commoditycode<=845019

replace isic4=27 if commoditycode>=940510 & commoditycode<=940560

replace isic4=27 if commoditycode>=850110 & commoditycode<=850790

replace isic4=27 if commoditycode>=850110 & commoditycode<=850790

replace isic4=27 if commoditycode>=850811 & commoditycode<=851390

replace isic4=27 if commoditycode>=851610 & commoditycode<=851690

replace isic4=27 if commoditycode>=853010 & commoditycode<=853190

replace isic4=27 if commoditycode>=853510 & commoditycode<=853990

replace isic4=27 if commoditycode>=854311 & commoditycode<=854590

replace isic4=27 if commoditycode>=854790 & commoditycode<=854890

********************(18)

replace isic4=28 if commoditycode==731511 | commoditycode==731512 | commoditycode==731519 | commoditycode==840130 | commoditycode==840120 | /// manufacture of machinery and equipment n.e.c
commoditycode==840721 | commoditycode==840729 | commoditycode==840790 | commoditycode==840810 | commoditycode==871620 | commoditycode==844391 | ///
commoditycode==840890 | commoditycode==841181 | commoditycode==841182 | commoditycode==841199 | commoditycode==848790 | ///
commoditycode==841459 | commoditycode==844311 | commoditycode==870110 | commoditycode==870130 | commoditycode==870190 | commoditycode==854310  ///

replace isic4=28 if commoditycode>=840510 & commoditycode<=840690

replace isic4=28 if commoditycode>=841011 & commoditycode<=841090

replace isic4=28 if commoditycode>=841480 & commoditycode<=841790

replace isic4=28 if commoditycode>=841221 & commoditycode<=841440

replace isic4=28 if commoditycode>=841850 & commoditycode<=841899

replace isic4=28 if commoditycode>=841931 & commoditycode<=842199

replace isic4=28 if commoditycode>=842220 & commoditycode<=844240

replace isic4=28 if commoditycode>=844313 & commoditycode<=844319

replace isic4=28 if commoditycode>=844400 & commoditycode<=844900

replace isic4=28 if commoditycode>=845020 & commoditycode<=845110

replace isic4=28 if commoditycode>=845129 & commoditycode<=846890

replace isic4=28 if commoditycode>=847410 & commoditycode<=848690

replace isic4=28 if commoditycode>=850810 & commoditycode<=850890

replace isic4=28 if commoditycode>=851410 & commoditycode<=851590

replace isic4=28 if commoditycode>=870911 & commoditycode<=871000

replace isic4=28 if twodigit==93

******************(19)

replace isic4=29 if commoditycode==840820 | commoditycode==840991 | commoditycode==840999 | /// manufacture of motor vehicles, trailers and semi-trailers
commoditycode==870120 | commoditycode==871610 | commoditycode==871690

replace isic4=29 if commoditycode>=840731 & commoditycode<=840734

replace isic4=29 if commoditycode>=870210 & commoditycode<=870899

replace isic4=29 if commoditycode>=871631 & commoditycode<=871640

*****************(20)

replace isic4=30 if commoditycode==840710 | commoditycode==871680 | commoditycode==840910 | /// manufacture of other transport equipment
commoditycode==841191 | commoditycode==841210

replace isic4=30 if twodigit==86

replace isic4=30 if twodigit==89

replace isic4=30 if twodigit==88 & commoditycode!=880400

replace isic4=30 if commoditycode>=841111 & commoditycode<=841122

replace isic4=30 if commoditycode>=871110 & commoditycode<=871499

****************(21)

replace isic4=31 if commoditycode>=940110 & commoditycode<=940190 /// manufacture of furniture

replace isic4=31 if commoditycode>=940310 & commoditycode<=940429

****************(22)

replace isic4=32 if twodigit==66 | twodigit==67 |  twodigit==92 | twodigit==95 /// other manufacturing 

replace isic4=32 if twodigit==96 & twodigit!=15

replace isic4=32 if commoditycode==340600 | commoditycode==360500 | commoditycode==360610 | commoditycode==360690 | commoditycode==420600 | ///
commoditycode==420610 | commoditycode==420690 | commoditycode==590410 | commoditycode==590491 | commoditycode==590492 | ///
commoditycode==420321 | commoditycode==871500 | commoditycode==902300 | commoditycode==711711 | commoditycode==711719 | commoditycode==711790

replace isic4=32 if commoditycode>=710122 & commoditycode<=710399

replace isic4=32 if commoditycode>=710490 & commoditycode<=710590

replace isic4=32 if commoditycode>=711311 & commoditycode<=711890

drop if commoditycode>=60110 & commoditycode<=60499 /// agriculture

drop if commoditycode>=130110 & commoditycode<=130239 /// forestry 

drop if commoditycode>=270111 & commoditycode<=270300 /// mining

drop if commoditycode==270500 | commoditycode==270900 | commoditycode==271410 | commoditycode==271490 | commoditycode==411520

drop if commoditycode>=391510 & commoditycode<=391590 /// goods not specified elsewhere

drop if commoditycode>=711230 & commoditycode<=711299 /// goods not specified elsewhere

drop if commoditycode>=400110 & commoditycode<=400130 /// agriculture

drop if commoditycode>=410320 & commoditycode<=410390 /// farming

drop if commoditycode==710110 | commoditycode==710121 

drop if twodigit==97 | twodigit==99

codebook isic4

**Only for 2015 and 2016 

drop if commoditycode==271600

*drop if isic4==.

************************************************************************************************************
