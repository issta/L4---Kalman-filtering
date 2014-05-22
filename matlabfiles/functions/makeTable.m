function makeTable( firstColumn,secondColumn,thirdColumn,fourthColumn,fifthColumn,sixthColumn,seventhColumn,eightColumn,tableName,titleTable,captionThing )
%makeTableLatex( firstColumn,secondColumn,thirdColumn,fourthColumn,tableName )
% Summary of this function goes here
%   Detailed explanation goes here
    currentD = cd;
	cd('/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L4 - Kalman filtering/Tables/')
    FID = fopen(tableName, 'w');
    fprintf(FID, '\\begin{table}[h] \n');
    fprintf(FID, '   \\begin{center} \n');
    fprintf(FID, '      \\begin{tabular}{lllllllll}\\toprule \n');
    fprintf(FID, ['\\multicolumn{5}{c}{' titleTable '}\\\\ \n']);
    fprintf(FID,'\\midrule \n');
    fprintf(FID,[...
        't',...
        '&  $x_e \\pm \\sigma_e$',...
        '& $x_n \\pm \\sigma_n$',...
        '& $v_e \\pm \\sigma_{v_e}$',...
        '& $v_n \\pm \\sigma_{v_n}$',...
        '\\\\ \\midrule \n']);
    t = 0:2:2*length(firstColumn);
    for k=1:length(firstColumn)
        fprintf(FID, '%g & $%8.2f \\pm %8.2f$ & $%8.2f \\pm %8.2f$ & $%8.2f \\pm %8.2f$ & $%8.2f \\pm %8.2f$  \\\\ ',t(k), firstColumn(k), secondColumn(k), thirdColumn(k), fourthColumn(k), fifthColumn(k),sixthColumn(k),seventhColumn(k),eightColumn(k));
        if k==length(firstColumn)
            fprintf(FID, '\\bottomrule ');
        end
        fprintf(FID, '\n');
    end
    fprintf(FID, '      \\end{tabular} \n');
    fprintf(FID, '   \\end{center}\n');
    fprintf(FID,['\\caption{' captionThing '} \n']);
    fprintf(FID,['\\label{tab:' tableName '} \n']);
    fprintf(FID, '\\end{table} \n');
    fclose(FID);
cd(currentD)
end

