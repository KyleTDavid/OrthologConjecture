import sys
from ete3 import Tree

#first argument is masterdata file, IN=paralog with a single species

tree=Tree("(((((((((((Poecilia_formosa:0,Xiphophorus_maculatus:0)A:4.070287,Oryzias_latipes:0)B:3.358582,Gasterosteus_aculeatus:0)C:3.05321,Oreochromis_niloticus:0)D:3.17584,(Takifugu_rubripes:0,Tetraodon_nigroviridis:0)E:3.968763)F:3.12639,Gadus_morhua:0)G:2.63493,(Astyanax_mexicanus:0,Danio_rerio:0)H:2.57855)I:1.76027,Lepisosteus_oculatus:0)J:0.81233,((((((((Dasypus_novemcinctus:0,Choloepus_hoffmanni:0)K:3.634505,(Echinops_telfairi:0,(Loxodonta_africana:0,Procavia_capensis:0)L:3.725591)M:3.621928)N:3.457881,(((Erinaceus_europaeus:0,Sorex_araneus:0)O:3.71388,((Felis_catus:0,(Canis_familiaris:0,(Ailuropoda_melanoleuca:0,Mustela_putorius_furo:0)P:3.873823)Q:3.868816)R:3.786047,((Equus_caballus:0,(Pteropus_vampyrus:0,Myotis_lucifugus:0)S:3.786597)T:3.611789,(Vicugna_pacos:0,(Sus_scrofa:0,(Tursiops_truncatus:0,(Bos_taurus:0,Ovis_aries:0)U:4.146878)V:3.800792)W:3.769051)X:3.71203)Y:3.598957)Z:3.601635)AA:3.561725,(((Oryctolagus_cuniculus:0,Ochotona_princeps:0)AB:3.848933,(Ictidomys_tridecemlineatus:0,((Dipodomys_ordii:0,(Jaculus_jaculus:0,(Nannospalax_galili:0,((Rattus_norvegicus:0,(Mus_pahari:0,(Mus_caroli:0,(Mus_spretus_spreteij:0,Mus_musculus:0)AC:4.2952388)AD:4.2696964)AE:4.2368343)AF:4.186918,(Peromyscus_maniculatus_bairdii:0,(Microtus_ochrogaster:0,(Mesocricetus_auratus:0,Cricetulus_griseus:0)AG:4.169692)AH:4.123307)AI:4.103811)AJ:4.052888)AK:3.968213)AL:3.863083)AM:3.730709,((Fukomys_damarensis:0,Heterocephalus_glaber:0)AN:4.096977,((Chinchilla_lanigera:0,Octodon_degus:0)AO:4.072171,(Cavia_porcellus:0,Cavia_aperea:0)AP:4.194472)AQ:4.041556)AR:3.987868)AS:3.677059)AT:3.686123)AU:3.618969,(((Microcebus_murinus:0,Otolemur_garnettii:0)AW:3.777662,(Carlito_syrichta:0,((Callithrix_jacchus:0),(((((Pan_troglodytes:0,Homo_sapiens:0)AY:4.2458682,Gorilla_gorilla:0)AZ:4.221375,Pongo_abelii:0)BA:4.153728,Nomascus_leucogenys:0)BB:4.12175,(Chlorocebus_sabaeus:0,(Macaca_mulatta:0,Papio_anubis:0)BC:4.222617)BD:4.188427)BE:4.053516)BF:3.939933)BG:3.676308)BH:3.645774,Tupaia_belangeri:0)BI:3.589621)BJ:0.81233)BK:3.469509)BL:3.32692,((Sarcophilus_harrisii:0,Notamacropus_eugenii:0)BM:3.816754,Monodelphis_domestica:0)BN:3.530633)BO:2.59139,Ornithorhynchus_anatinus:0)BP:2.40478,(((((Gallus_gallus:0,Meleagris_gallopavo:0)BQ:4.175152,Anas_platyrhynchos:0)BR:3.781495,(Ficedula_albicollis:0,Taeniopygia_guttata:0)BS:4.088268)BT:3.59737,Pelodiscus_sinensis:0)BU:2.07443,Anolis_carolinensis:0)BV:1.22813)BW:1.37177,Xenopus_tropicalis:0)BX:0.83961,Latimeria_chalumnae:0)BY:0.18859)AV):1,Petromyzon_marinus:0)XX;",format=1)

with open(sys.argv[1]) as f:
    for line in f:
        done = []
        #get genes in each node
        linelist=line.split('\t')
        linelistA = linelist[9].split(', ')
        linelistB = linelist[14].split(', ')
        linelistC = linelistA + linelistB
        for n, i in enumerate(linelistC):
            linelistC[n] = ''.join([x for x in i if not x.isdigit()])
            if linelistC[n] == "ENSACAP":
                linelistC[n] = "Anolis_carolinensis"
            if linelistC[n] == "ENSAMEP":
                linelistC[n] = "Ailuropoda_melanoleuca"
            if linelistC[n] == "ENSAMXP":
                linelistC[n] = "Astyanax_mexicanus"
            if linelistC[n] == "ENSAPLP":
                linelistC[n] = "Anas_platyrhynchos"
            if linelistC[n] == "ENSBTAP":
                linelistC[n] = "Bos_taurus"
            if linelistC[n] == "ENSCAFP":
                linelistC[n] = "Canis_familiaris"
            if linelistC[n] == "ENSCAPP":
                linelistC[n] = "Cavia_aperea"
            if linelistC[n] == "ENSCGRP":
                linelistC[n] = "Cricetulus_griseus"
            if linelistC[n] == "ENSCHOP":
                linelistC[n] = "Choloepus_hoffmanni"
            if linelistC[n] == "ENSCJAP":
                linelistC[n] = "Callithrix_jacchus"
            if linelistC[n] == "ENSCLAP":
                linelistC[n] = "Chinchilla_lanigera"
            if linelistC[n] == "ENSCPOP":
                linelistC[n] = "Cavia_porcellus"
            if linelistC[n] == "ENSCSAP":
                linelistC[n] = "Chlorocebus_sabaeus"
            if linelistC[n] == "ENSDARP":
                linelistC[n] = "Danio_rerio"
            if linelistC[n] == "ENSDNOP":
                linelistC[n] = "Dasypus_novemcinctus"
            if linelistC[n] == "ENSDORP":
                linelistC[n] = "Dipodomys_ordii"
            if linelistC[n] == "ENSECAP":
                linelistC[n] = "Equus_caballus"
            if linelistC[n] == "ENSEEUP":
                linelistC[n] = "Erinaceus_europaeus"
            if linelistC[n] == "ENSETEP":
                linelistC[n] = "Echinops_telfairi"
            if linelistC[n] == "ENSFALP":
                linelistC[n] = "Ficedula_albicollis"
            if linelistC[n] == "ENSFCAP":
                linelistC[n] = "Felis_catus"
            if linelistC[n] == "ENSFDAP":
                linelistC[n] = "Fukomys_damarensis"
            if linelistC[n] == "ENSGACP":
                linelistC[n] = "Gasterosteus_aculeatus"
            if linelistC[n] == "ENSGALP":
                linelistC[n] = "Gallus_gallus"
            if linelistC[n] == "ENSGGOP":
                linelistC[n] = "Gorilla_gorilla"
            if linelistC[n] == "ENSGMOP":
                linelistC[n] = "Gadus_morhua"
            if linelistC[n] == "ENSHGLP":
                linelistC[n] = "Heterocephalus_glaber"
            if linelistC[n] == "ENSJJAP":
                linelistC[n] = "Jaculus_jaculus"
            if linelistC[n] == "ENSLACP":
                linelistC[n] = "Latimeria_chalumnae"
            if linelistC[n] == "ENSLAFP":
                linelistC[n] = "Loxodonta_africana"
            if linelistC[n] == "ENSLOCP":
                linelistC[n] = "Lepisosteus_oculatus"
            if linelistC[n] == "ENSMAUP":
                linelistC[n] = "Mesocricetus_auratus"
            if linelistC[n] == "ENSMEUP":
                linelistC[n] = "Notamacropus_eugenii"
            if linelistC[n] == "ENSMGAP":
                linelistC[n] = "Meleagris_gallopavo"
            if linelistC[n] == "ENSMICP":
                linelistC[n] = "Microcebus_murinus"
            if linelistC[n] == "ENSMLUP":
                linelistC[n] = "Myotis_lucifugus"
            if linelistC[n] == "ENSMMUP":
                linelistC[n] = "Macaca_mulatta"
            if linelistC[n] == "ENSMOCP":
                linelistC[n] = "Microtus_ochrogaster"
            if linelistC[n] == "ENSMODP":
                linelistC[n] = "Monodelphis_domestica"
            if linelistC[n] == "ENSMPUP":
                linelistC[n] = "Mustela_putorius_furo"
            if linelistC[n] == "ENSMUSP":
                linelistC[n] = "Mus_musculus"
            if linelistC[n] == "ENSNGAP":
                linelistC[n] = "Nannospalax_galili"
            if linelistC[n] == "ENSNLEP":
                linelistC[n] = "Nomascus_leucogenys"
            if linelistC[n] == "ENSOANP":
                linelistC[n] = "Ornithorhynchus_anatinus"
            if linelistC[n] == "ENSOARP":
                linelistC[n] = "Ovis_aries"
            if linelistC[n] == "ENSOCUP":
                linelistC[n] = "Oryctolagus_cuniculus"
            if linelistC[n] == "ENSODEP":
                linelistC[n] = "Octodon_degus"
            if linelistC[n] == "ENSOGAP":
                linelistC[n] = "Otolemur_garnettii"
            if linelistC[n] == "ENSONIP":
                linelistC[n] = "Oreochromis_niloticus"
            if linelistC[n] == "ENSOPRP":
                linelistC[n] = "Ochotona_princeps"
            if linelistC[n] == "ENSORLP":
                linelistC[n] = "Oryzias_latipes"
            if linelistC[n] == "ENSP":
                linelistC[n] = "Homo_sapiens"
            if linelistC[n] == "ENSPANP":
                linelistC[n] = "Papio_anubis"
            if linelistC[n] == "ENSPCAP":
                linelistC[n] = "Procavia_capensis"
            if linelistC[n] == "ENSPEMP":
                linelistC[n] = "Peromyscus_maniculatus_bairdii"
            if linelistC[n] == "ENSPFOP":
                linelistC[n] = "Poecilia_formosa"
            if linelistC[n] == "ENSPPYP":
                linelistC[n] = "Pongo_abelii"
            if linelistC[n] == "ENSPSIP":
                linelistC[n] = "Pelodiscus_sinensis"
            if linelistC[n] == "ENSPTRP":
                linelistC[n] = "Pan_troglodytes"
            if linelistC[n] == "ENSPVAP":
                linelistC[n] = "Pteropus_vampyrus"
            if linelistC[n] == "ENSRNOP":
                linelistC[n] = "Rattus_norvegicus"
            if linelistC[n] == "ENSSARP":
                linelistC[n] = "Sorex_araneus"
            if linelistC[n] == "ENSSHAP":
                linelistC[n] = "Sarcophilus_harrisii"
            if linelistC[n] == "ENSSSCP":
                linelistC[n] = "Sus_scrofa"
            if linelistC[n] == "ENSSTOP":
                linelistC[n] = "Ictidomys_tridecemlineatus"
            if linelistC[n] == "ENSTBEP":
                linelistC[n] = "Tupaia_belangeri"
            if linelistC[n] == "ENSTGUP":
                linelistC[n] = "Taeniopygia_guttata"
            if linelistC[n] == "ENSTNIP":
                linelistC[n] = "Tetraodon_nigroviridis"
            if linelistC[n] == "ENSTRUP":
                linelistC[n] = "Takifugu_rubripes"
            if linelistC[n] == "ENSTSYP":
                linelistC[n] = "Carlito_syrichta"
            if linelistC[n] == "ENSTTRP":
                linelistC[n] = "Tursiops_truncatus"
            if linelistC[n] == "ENSVPAP":
                linelistC[n] = "Vicugna_pacos"
            if linelistC[n] == "ENSXETP":
                linelistC[n] = "Xenopus_tropicalis"
            if linelistC[n] == "ENSXMAP":
                linelistC[n] = "Xiphophorus_maculatus"
            if linelistC[n] == "MGP_CAROLIEiJ_P":
                linelistC[n] = "Mus_caroli"
            if linelistC[n] == "MGP_PahariEiJ_P":
                linelistC[n] = "Mus_pahari"
            if linelistC[n] == "MGP_SPRETEiJ_P":
                linelistC[n] = "Mus_spretus_spreteij"
            if linelistC[n] == "ENSPMAP":
                linelistC[n] = "Petromyzon_marinus"

        if len(linelistC)==2 and linelistC[0] == linelistC[1]:
            print "IN"
        else:
            MRCA=tree.get_common_ancestor(linelistC)
            if MRCA.is_leaf() == False:
                print MRCA.name






