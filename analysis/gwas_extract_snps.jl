Distributed.@everywhere using RCall,CSV,DataFrames,GZip,Random
Distributed.@everywhere include("gwas_tools.jl")

function call_gwas_extract_snps(somaid_arr,fdr_level,root_path,on_google,study)

    data_folder = "$(root_path)/data/gwas_protein/"
    out_folder = "$(root_path)/data/gwas_protein_out/"

    if study == "jhs"
        #JHS
        #gwas_file_pre = "jhs_prot_adjust_age_sex_pcs_prot_"
        gwas_file_pre = ""
        gwas_file_post = ".fastgwa.gz"
    elseif study == "mesa"
        #MESA
        gwas_file_pre = ""
        gwas_file_post = "_1_resid.txt_baseline.fastgwa.gz"
    elseif study == "heritage"
        #HERITAGE
        gwas_file_pre = "B_"
        gwas_file_post = ".txt.fastgwa.gz"
    else
        @error("Study $(study) not recognized!")
    end

    gwas_file_keep_post = ".fastgwa.txt"

    #data_folder = "/Users/aaroneisman/projects/benson/data/gwas_protein_out/"
    #somaid_arr = ["SL005699","SL012707","SL000020"]

    #gwas_fh = "/Users/aaroneisman/projects/benson/data/dan/jhs_prot_adjust_age_sex_pcs_prot_SL005699.fastgwa.gz"
    #gwas_keep_fh = "gwas_keep_snps.txt"

    #fdr_level = 0.1
    #keep_file = out_folder*gwas_file_pre*"all_fdr_$(fdr_level).txt"
    keep_file = out_folder*"/ldr/all_ldr_$(fdr_level).txt"
    keep_df = CSV.read(keep_file,DataFrame)
    keep_snp_set = Set(keep_df.SNP)

    shuffle!(somaid_arr)

    Distributed.@distributed vcat for somaid in somaid_arr
        if on_google == 1
            ## EXAMPLE: run(`gsutil cp gs://jhs_data_topmed/jhs_prot_adjust_age_sex_pcs/prot/SL000001.fastgwa.gz SL000001.fastgwa.gz​`)
            gwas_file = "$(gwas_file_pre)$(somaid)$(gwas_file_post)"
            sleep(rand()*10)
            if study == "jhs"
                #JHS
                run(`gsutil cp gs://jhs_data_topmed/jhs_gwas/jhs_SOMA_gwas/fastgwa_with_pcs/raw_results/$(gwas_file) $(data_folder)$(gwas_file)`)
            elseif study == "mesa"
                #MESA
                run(`gsutil cp gs://mesa-bucket/mesa_gwas/mesa_SOMA_gwas/fastgwa/raw_results/$(gwas_file) $(data_folder)$(gwas_file)`)
            elseif study == "heritage"
                #HERITAGE
                run(`gsutil cp gs://heritage-bucket/heritage_gwas/heritage_SOMA_gwas/fastgwa/raw_results/$(gwas_file) $(data_folder)$(gwas_file)`)
            else
                @error("Study $(study) not recognized")
            end
        end
        
        try
            gwas_fh = data_folder*gwas_file_pre*somaid*gwas_file_post
            gwas_keep_fh = out_folder*gwas_file_pre*somaid*"_keep"*gwas_file_keep_post
            gwas_extract_snps(gwas_fh,gwas_keep_fh,keep_snp_set,"\t")
        catch
            @warn("Could not extract snps for $(gwas_fh)")
        end
        
        if on_google == 1
            ## EXAMPLE: run(`gsutil cp gs://jhs_data_topmed/jhs_prot_adjust_age_sex_pcs/prot/SL000001.fastgwa.gz SL000001.fastgwa.gz​`)
            gwas_file = "$(gwas_file_pre)$(somaid)$(gwas_file_post)"
            run(`rm $(data_folder)$(gwas_file)`)
        end
    end

end

#root_path = "/home/aaron_eisman"
#root_path = "/home/aaroneisman"
#on_google = 1
#fdr_level = 0.05

#somaid_arr = ["SL000001","SL000002","SL000003","SL000004","SL000006","SL000007","SL000009","SL000011","SL000017","SL000019","SL000020","SL000021","SL000022","SL000024","SL000027","SL000038","SL000039","SL000045","SL000047","SL000048","SL000049","SL000051","SL000052","SL000053","SL000055","SL000057","SL000062","SL000064","SL000070","SL000074","SL000076","SL000084","SL000087","SL000088","SL000089","SL000104","SL000124","SL000125","SL000130","SL000131","SL000133","SL000134","SL000136","SL000137","SL000138","SL000139","SL000142","SL000145","SL000158","SL000164","SL000247","SL000248","SL000249","SL000250","SL000251","SL000252","SL000254","SL000268","SL000271","SL000272","SL000276","SL000277","SL000280","SL000283","SL000299","SL000300","SL000305","SL000306","SL000308","SL000309","SL000310","SL000312","SL000313","SL000314","SL000316","SL000318","SL000319","SL000320","SL000321","SL000322","SL000323","SL000324","SL000325","SL000337","SL000338","SL000339","SL000342","SL000343","SL000344","SL000345","SL000346","SL000347","SL000357","SL000358","SL000360","SL000377","SL000382","SL000383","SL000384","SL000396","SL000398","SL000401","SL000403","SL000406","SL000408","SL000409","SL000414","SL000415","SL000420","SL000424","SL000426","SL000427","SL000428","SL000433","SL000437","SL000441","SL000442","SL000445","SL000448","SL000449","SL000450","SL000451","SL000454","SL000455","SL000456","SL000458","SL000459","SL000460","SL000461","SL000462","SL000466","SL000467","SL000468","SL000470","SL000474","SL000478","SL000479","SL000480","SL000481","SL000483","SL000485","SL000493","SL000496","SL000497","SL000498","SL000506","SL000507","SL000508","SL000509","SL000510","SL000515","SL000516","SL000517","SL000519","SL000521","SL000522","SL000523","SL000524","SL000525","SL000526","SL000527","SL000528","SL000530","SL000532","SL000535","SL000537","SL000539","SL000540","SL000541","SL000542","SL000545","SL000546","SL000550","SL000551","SL000553","SL000554","SL000556","SL000557","SL000558","SL000560","SL000563","SL000565","SL000566","SL000570","SL000572","SL000573","SL000581","SL000582","SL000584","SL000586","SL000587","SL000588","SL000589","SL000590","SL000591","SL000592","SL000597","SL000598","SL000601","SL000603","SL000605","SL000613","SL000615","SL000616","SL000617","SL000622","SL000633","SL000638","SL000640","SL000645","SL000655","SL000658","SL000668","SL000670","SL000674","SL000678","SL000695","SL000836","SL001691","SL001713","SL001716","SL001717","SL001718","SL001720","SL001721","SL001726","SL001728","SL001729","SL001731","SL001737","SL001753","SL001761","SL001766","SL001774","SL001777","SL001791","SL001795","SL001796","SL001797","SL001800","SL001802","SL001815","SL001880","SL001888","SL001890","SL001896","SL001897","SL001902","SL001905","SL001938","SL001943","SL001945","SL001947","SL001973","SL001990","SL001992","SL001995","SL001996","SL001997","SL001998","SL001999","SL002036","SL002071","SL002075","SL002078","SL002081","SL002086","SL002093","SL002505","SL002506","SL002508","SL002517","SL002519","SL002522","SL002524","SL002525","SL002528","SL002539","SL002541","SL002542","SL002561","SL002565","SL002602","SL002621","SL002640","SL002644","SL002646","SL002650","SL002654","SL002655","SL002662","SL002684","SL002688","SL002695","SL002702","SL002704","SL002705","SL002706","SL002722","SL002731","SL002755","SL002756","SL002762","SL002763","SL002782","SL002783","SL002785","SL002792","SL002803","SL002823","SL002922","SL003041","SL003043","SL003060","SL003066","SL003080","SL003104","SL003166","SL003167","SL003168","SL003169","SL003170","SL003171","SL003172","SL003173","SL003176","SL003177","SL003178","SL003179","SL003182","SL003183","SL003184","SL003186","SL003187","SL003188","SL003189","SL003190","SL003191","SL003192","SL003193","SL003196","SL003197","SL003198","SL003199","SL003200","SL003201","SL003220","SL003280","SL003300","SL003301","SL003302","SL003303","SL003304","SL003305","SL003307","SL003308","SL003309","SL003310","SL003320","SL003322","SL003323","SL003324","SL003326","SL003327","SL003328","SL003329","SL003331","SL003332","SL003334","SL003340","SL003341","SL003349","SL003362","SL003440","SL003461","SL003520","SL003522","SL003524","SL003542","SL003643","SL003646","SL003647","SL003648","SL003650","SL003653","SL003655","SL003657","SL003658","SL003672","SL003674","SL003679","SL003680","SL003685","SL003687","SL003690","SL003700","SL003703","SL003704","SL003710","SL003711","SL003717","SL003722","SL003726","SL003728","SL003733","SL003735","SL003738","SL003739","SL003744","SL003753","SL003755","SL003761","SL003764","SL003770","SL003774","SL003785","SL003792","SL003793","SL003800","SL003803","SL003848","SL003849","SL003862","SL003863","SL003869","SL003872","SL003915","SL003916","SL003918","SL003919","SL003930","SL003951","SL003970","SL003974","SL003990","SL003993","SL003994","SL004008","SL004009","SL004010","SL004015","SL004016","SL004060","SL004063","SL004064","SL004066","SL004067","SL004068","SL004070","SL004078","SL004080","SL004081","SL004097","SL004101","SL004118","SL004119","SL004120","SL004125","SL004126","SL004128","SL004131","SL004133","SL004134","SL004136","SL004137","SL004139","SL004140","SL004141","SL004142","SL004143","SL004144","SL004145","SL004146","SL004147","SL004148","SL004149","SL004151","SL004152","SL004153","SL004154","SL004155","SL004156","SL004157","SL004159","SL004160","SL004180","SL004182","SL004183","SL004208","SL004209","SL004230","SL004248","SL004253","SL004258","SL004260","SL004269","SL004271","SL004296","SL004297","SL004298","SL004299","SL004301","SL004304","SL004305","SL004306","SL004326","SL004327","SL004329","SL004330","SL004331","SL004332","SL004333","SL004334","SL004335","SL004336","SL004337","SL004338","SL004339","SL004340","SL004342","SL004343","SL004345","SL004346","SL004347","SL004348","SL004349","SL004350","SL004351","SL004352","SL004353","SL004354","SL004355","SL004356","SL004357","SL004359","SL004360","SL004362","SL004363","SL004364","SL004365","SL004366","SL004367","SL004396","SL004400","SL004415","SL004438","SL004457","SL004458","SL004466","SL004467","SL004469","SL004475","SL004477","SL004482","SL004484","SL004486","SL004489","SL004492","SL004511","SL004515","SL004516","SL004536","SL004556","SL004557","SL004579","SL004580","SL004588","SL004589","SL004591","SL004594","SL004605","SL004610","SL004625","SL004626","SL004635","SL004636","SL004637","SL004639","SL004642","SL004643","SL004644","SL004645","SL004646","SL004648","SL004649","SL004650","SL004652","SL004654","SL004660","SL004661","SL004668","SL004669","SL004670","SL004671","SL004672","SL004673","SL004676","SL004683","SL004685","SL004686","SL004687","SL004689","SL004690","SL004692","SL004697","SL004698","SL004704","SL004708","SL004712","SL004714","SL004716","SL004718","SL004720","SL004723","SL004724","SL004725","SL004726","SL004733","SL004737","SL004739","SL004742","SL004747","SL004750","SL004751","SL004752","SL004757","SL004759","SL004760","SL004765","SL004768","SL004771","SL004781","SL004782","SL004783","SL004795","SL004804","SL004805","SL004811","SL004812","SL004814","SL004815","SL004820","SL004821","SL004823","SL004827","SL004829","SL004837","SL004838","SL004843","SL004844","SL004845","SL004849","SL004850","SL004851","SL004852","SL004853","SL004855","SL004856","SL004858","SL004859","SL004860","SL004861","SL004862","SL004863","SL004864","SL004865","SL004866","SL004867","SL004868","SL004869","SL004871","SL004872","SL004875","SL004876","SL004891","SL004899","SL004901","SL004908","SL004910","SL004914","SL004915","SL004919","SL004920","SL004921","SL004924","SL004925","SL004932","SL004938","SL004939","SL004940","SL004984","SL005034","SL005059","SL005084","SL005087","SL005102","SL005115","SL005152","SL005153","SL005155","SL005156","SL005157","SL005158","SL005159","SL005160","SL005161","SL005164","SL005165","SL005166","SL005167","SL005168","SL005169","SL005170","SL005171","SL005172","SL005173","SL005174","SL005177","SL005178","SL005179","SL005181","SL005183","SL005184","SL005185","SL005187","SL005188","SL005189","SL005190","SL005191","SL005193","SL005194","SL005195","SL005196","SL005197","SL005199","SL005200","SL005201","SL005202","SL005204","SL005205","SL005206","SL005207","SL005208","SL005209","SL005210","SL005212","SL005213","SL005214","SL005215","SL005217","SL005218","SL005219","SL005220","SL005221","SL005222","SL005223","SL005224","SL005225","SL005226","SL005227","SL005228","SL005229","SL005230","SL005231","SL005233","SL005234","SL005235","SL005236","SL005250","SL005256","SL005258","SL005261","SL005263","SL005266","SL005308","SL005352","SL005357","SL005358","SL005361","SL005372","SL005392","SL005403","SL005430","SL005437","SL005488","SL005491","SL005493","SL005508","SL005572","SL005574","SL005575","SL005588","SL005629","SL005630","SL005675","SL005679","SL005685","SL005687","SL005688","SL005694","SL005699","SL005703","SL005725","SL005730","SL005764","SL005789","SL005793","SL005797","SL005846","SL006029","SL006088","SL006091","SL006108","SL006114","SL006119","SL006131","SL006132","SL006189","SL006197","SL006230","SL006268","SL006372","SL006374","SL006378","SL006397","SL006406","SL006448","SL006460","SL006476","SL006480","SL006512","SL006522","SL006523","SL006527","SL006528","SL006542","SL006544","SL006550","SL006610","SL006629","SL006675","SL006694","SL006698","SL006705","SL006713","SL006777","SL006803","SL006805","SL006830","SL006892","SL006910","SL006911","SL006912","SL006913","SL006914","SL006915","SL006916","SL006917","SL006918","SL006919","SL006920","SL006921","SL006922","SL006923","SL006924","SL006970","SL006992","SL006993","SL006998","SL007003","SL007022","SL007024","SL007025","SL007033","SL007049","SL007056","SL007059","SL007070","SL007100","SL007108","SL007121","SL007122","SL007136","SL007145","SL007151","SL007153","SL007173","SL007179","SL007181","SL007195","SL007206","SL007207","SL007217","SL007221","SL007223","SL007228","SL007229","SL007237","SL007242","SL007250","SL007261","SL007266","SL007272","SL007274","SL007280","SL007281","SL007284","SL007295","SL007306","SL007310","SL007311","SL007324","SL007327","SL007328","SL007336","SL007356","SL007358","SL007361","SL007373","SL007385","SL007429","SL007453","SL007471","SL007502","SL007531","SL007547","SL007560","SL007620","SL007631","SL007640","SL007642","SL007651","SL007673","SL007674","SL007680","SL007696","SL007729","SL007747","SL007752","SL007756","SL007774","SL007804","SL007806","SL007828","SL007869","SL007871","SL007888","SL007953","SL008008","SL008023","SL008039","SL008059","SL008063","SL008071","SL008072","SL008085","SL008094","SL008099","SL008102","SL008113","SL008122","SL008143","SL008157","SL008158","SL008176","SL008177","SL008178","SL008190","SL008193","SL008331","SL008360","SL008372","SL008378","SL008380","SL008381","SL008382","SL008402","SL008414","SL008416","SL008421","SL008466","SL008486","SL008504","SL008516","SL008522","SL008574","SL008588","SL008590","SL008591","SL008609","SL008611","SL008614","SL008623","SL008631","SL008639","SL008703","SL008709","SL008728","SL008759","SL008760","SL008773","SL008808","SL008810","SL008822","SL008835","SL008837","SL008865","SL008904","SL008909","SL008916","SL008931","SL008933","SL008936","SL008945","SL008956","SL008967","SL009045","SL009054","SL009089","SL009202","SL009207","SL009210","SL009213","SL009216","SL009324","SL009328","SL009341","SL009349","SL009400","SL009412","SL009431","SL009628","SL009629","SL009768","SL009790","SL009791","SL009792","SL009868","SL009948","SL009951","SL009988","SL010250","SL010288","SL010328","SL010348","SL010349","SL010368","SL010369","SL010371","SL010372","SL010373","SL010374","SL010375","SL010376","SL010378","SL010379","SL010381","SL010384","SL010388","SL010390","SL010391","SL010393","SL010449","SL010450","SL010451","SL010454","SL010455","SL010456","SL010457","SL010458","SL010459","SL010460","SL010461","SL010462","SL010463","SL010464","SL010465","SL010466","SL010467","SL010468","SL010469","SL010470","SL010471","SL010488","SL010489","SL010490","SL010491","SL010492","SL010493","SL010494","SL010495","SL010496","SL010498","SL010499","SL010500","SL010501","SL010502","SL010503","SL010504","SL010505","SL010508","SL010509","SL010510","SL010512","SL010513","SL010514","SL010515","SL010516","SL010517","SL010518","SL010519","SL010520","SL010521","SL010522","SL010523","SL010524","SL010528","SL010529","SL010530","SL010610","SL010612","SL010613","SL010616","SL010617","SL010619","SL010830","SL010927","SL010928","SL010973","SL011049","SL011068","SL011069","SL011071","SL011073","SL011100","SL011180","SL011202","SL011211","SL011232","SL011400","SL011404","SL011405","SL011406","SL011448","SL011498","SL011499","SL011508","SL011509","SL011510","SL011528","SL011529","SL011530","SL011532","SL011533","SL011535","SL011549","SL011588","SL011616","SL011628","SL011629","SL011630","SL011631","SL011708","SL011709","SL011768","SL011769","SL011770","SL011772","SL011808","SL011809","SL011888","SL012108","SL012148","SL012168","SL012188","SL012248","SL012395","SL012457","SL012469","SL012517","SL012538","SL012561","SL012698","SL012707","SL012740","SL012754","SL012769","SL012774","SL012783","SL012822","SL012881","SL013165","SL013240","SL013488","SL013489","SL013490","SL013548","SL013570","SL013754","SL013872","SL013928","SL013969","SL013988","SL013989","SL014008","SL014009","SL014028","SL014029","SL014048","SL014069","SL014070","SL014071","SL014088","SL014091","SL014092","SL014093","SL014094","SL014096","SL014108","SL014111","SL014113","SL014129","SL014130","SL014148","SL014188","SL014208","SL014209","SL014228","SL014229","SL014248","SL014268","SL014269","SL014270","SL014288","SL014289","SL014292","SL014294","SL014308","SL014468","SL014469","SL014470","SL014488","SL014684","SL014735","SL014896","SL014983","SL015046","SL015510","SL015728","SL016128","SL016129","SL016130","SL016148","SL016548","SL016549","SL016550","SL016551","SL016553","SL016554","SL016555","SL016557","SL016563","SL016566","SL016567","SL016828","SL016928","SL016969","SL017106","SL017128","SL017188","SL017189","SL017289","SL017290","SL017328","SL017424","SL017528","SL017529","SL017610","SL017611","SL017612","SL017614","SL018256","SL018509","SL018548","SL018587","SL018625","SL018629","SL018887","SL018891","SL018900","SL018921","SL018938","SL018946","SL018947","SL018971","SL019019","SL019096","SL019100","SL019978","SL019979","SL020171","SL020172","SL021043"]

#call_gwas_extract_snps(somaid_arr,fdr_level,root_path,on_google)

