folder = ARGS[1]
ext = ARGS[2]


function combine(folder, ext)

    open("$(folder)/$(ext).comb","w") do io

    files_arr = readdir(folder, sort = false);
    ext_files = files_arr[endswith.(files_arr,ext)];    

        for file in ext_files
            println(file)
            for line in eachline("$(folder)/$(file)")
                write(io,file)
                write(io,",")
                write(io,folder)
                write(io,",")
                write(io,line)
                write(io,"\n")
            end
        end
        
    end

end

combine(folder, ext)
