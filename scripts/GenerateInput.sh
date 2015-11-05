# this script generates the input for different temperatures and different densities
for T in 1*10^6 1*10^7 1*10^8 
    do
        #generate the folder
        mkdir -p T$T 
        #copy the required files
        cp gen_pa_ilkka.py T$T/
        cp ilkka.xml T$T/
        cd T$T
        #calculate beta
        beta=$(echo "scale=10; 1.0/($T * 3.16681*10^(-6))" | bc -l)
        #change the important files
        perl -i -spe "s/beta=\"(.*?)\"/beta=\"$beta\"/" ilkka.xml
        perl -i -spe "s/beta = (.*?) #/beta = $beta #/" gen_pa_ilkka.py
        echo "Start generating the pair potentials for the system (T=$T)"
        python -E gen_pa_ilkka.py
        echo "ended the generation of the pair potentials (T=$T)"
        for rho in 1000 100 10 1 0.1 0.01 0.001
            do
                #generate the subfolder
                mkdir -p rho$rho
                cd rho$rho
                #cp the files from above, because the temperature is fixed, only the breakup has to be recalculated
                cp ../*.h5 .
                cp ../*.txt .
                cp ../*.dat .
                cp ../INPUT .
                cp ../gen_pa_ilkka.py .
                cp ../ilkka.xml .
                #Calculate rho in terms of particles per angstroem^3
                cm_div_angstroem=100000000 
                mBe_div_g=$(echo "scale=10; 1.496508*10^(-23)"| bc -l)
                rho2=$(echo "scale=10; 1.0/($mBe_div_g*($cm_div_angstroem^3))" | bc -l)
                NBeAtoms=1
                NElectrons=$(echo "scale=10; 3*$NBeAtoms" | bc -l)
                #Calculate L
                L=$(echo "scale=100; e(l($NBeAtoms/$rho2)/3)" | bc -l)
                echo $L
                #Write the things to the xml file
                #perl -i -spe "s/output_prefix=\".*\"/output_prefix=\"WDM.T$T.rho$rho\"/" ilkka.xml
                perl -i -spe "s/L=\"(.*?)\"/L=\"$L\"/" ilkka.xml
                perl -i -spe "s/name=\"e\" n_part=\"(.*?)\"/name=\"e\" n_part=\"$NElectrons\"/" ilkka.xml
                perl -i -spe "s/name=\"Be\" n_part=\"(.*?)\"/name=\"e\" n_part=\"$NBeAtoms\"/" ilkka.xml
                #write the things to the python file
                perl -i -spe "s/L = (.*?) #/L = $L #/" gen_pa_ilkka.py
                #python gen_pa_ilkka.py
                cd ..
            done
        cd ..
    done
