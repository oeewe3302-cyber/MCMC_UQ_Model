function batches = load_excel_data(filepath)
% LOAD_EXCEL_DATA Load experimental data from Excel file with unit conversion
%
% Syntax:
%   batches = load_excel_data(filepath)
%
% Input:
%   filepath - Path to Excel file
%
% Output:
%   batches - Cell array of batch structures
%
% Each batch structure contains:
%   .name           - Sheet name
%   .t              - Time (h)
%   .x              - Viable cell density (10^6 cells/mL)
%   .via            - Viability (fraction)
%   .glc            - Glucose concentration (mM) - CONVERTED from g/L
%   .gln            - Glutamine concentration (mM) - CONVERTED from g/L
%   .lac            - Lactate concentration (mM) - CONVERTED from g/L
%   .nh4            - Ammonia concentration (mM) - CONVERTED from g/L
%   .vol            - Volume (mL) - CONVERTED from L
%   .feed_glc_vol   - Pure glucose feed volume (mL) - CONVERTED from L
%   .feed_vols      - Feed media volumes [n_time x n_feeds] (mL) - CONVERTED from L
%   .n_feeds        - Number of feed media (1-5)
%   .feed_glc       - Glucose conc in each feed media [n_feeds x 1] (mM) - CONVERTED
%   .feed_gln       - Glutamine conc in each feed media [n_feeds x 1] (mM) - CONVERTED
%   .glc_stock      - Glucose stock solution conc (mM) - CONVERTED
%
% Unit Conversion:
%   Excel data is in g/L, converted to mM using molecular weights:
%   - Glucose (GLC):   MW = 180.16 g/mol  → mM = (g/L) / 180.16 * 1000
%   - Glutamine (GLN): MW = 146.14 g/mol  → mM = (g/L) / 146.14 * 1000
%   - Lactate (LAC):   MW = 90.08 g/mol   → mM = (g/L) / 90.08 * 1000
%   - Ammonia (NH4):   MW = 18.04 g/mol   → mM = (g/L) / 18.04 * 1000

    % ================================================
    % Molecular weights (g/mol)
    % ================================================
    MW_GLC = 180.16;  % Glucose
    MW_GLN = 146.14;  % Glutamine
    MW_LAC = 90.08;   % Lactate
    MW_NH4 = 18.04;   % Ammonia (NH4+)
    
    % ================================================
    % 1. Read Main parameter sheet
    % ================================================
    try
       T_param = readtable('data.xlsx', 'Sheet', 'Main parameter', ...
                            'VariableNamingRule', 'preserve', ...
                            'Range', 'B2:C7'); 
        % Extract feed compositions (rows 2-6 for feeds 1-5, row 7 for glucose stock)
        param_data = table2array(T_param);  % Skip header row
        
        % Feed media 1-5 (rows 1-5 of param_data) - IN g/L, CONVERT TO mM
        feed_glc_all_gL = param_data(1:5, 1);  % Column 2: Glucose (g/L)
        feed_gln_all_gL = param_data(1:5, 2);  % Column 3: Glutamine (g/L)
        
        % Convert to mM
        feed_glc_all = feed_glc_all_gL / MW_GLC * 1000;  % mM
        feed_gln_all = feed_gln_all_gL / MW_GLN * 1000;  % mM
        
        % Glucose stock solution (row 6 of param_data = row 7 of sheet) - IN g/L
        glc_stock_gL = param_data(6, 1);
        glc_stock = glc_stock_gL / MW_GLC * 1000;  % Convert to mM
        
        fprintf('✅ Main parameter sheet loaded:\n');
        fprintf('   Feed media compositions:\n');
        for i = 1:5
            fprintf('     Feed %d: GLC=%.1f g/L (%.2f mM), GLN=%.1f g/L (%.2f mM)\n', ...
                    i, feed_glc_all_gL(i), feed_glc_all(i), ...
                    feed_gln_all_gL(i), feed_gln_all(i));
        end
        fprintf('   Glucose stock: %.1f g/L (%.2f mM)\n\n', glc_stock_gL, glc_stock);
        
    catch ME
        error('Failed to read Main parameter sheet: %s', ME.message);
    end
    
    % ================================================
    % 2. Read batch data sheets
    % ================================================
    [~, sheets] = xlsfinfo('data.xlsx');
    
    % Remove 'Main parameter' from sheets list
    sheets = sheets(~strcmp(sheets, 'Main parameter'));
    
    nB = numel(sheets);
    batches = cell(1, nB);

    for i = 1:nB
        fprintf('Loading batch %d/%d: %s...\n', i, nB, sheets{i});
        
        T = readtable('data.xlsx', 'Sheet', sheets{i}, 'VariableNamingRule', 'preserve');
        T = standardizeMissing(T, {'', 'NA', '#N/A', 'n/a', 'NaN', '-'});
        
        % Remove header row if exists
        if any(cellfun(@ischar, table2cell(T(1,:))))
            T(1,:) = []; 
        end
        
        raw = table2array(T);
        if iscell(raw)
            raw = cellfun(@toNum, raw); 
        end
        
        % Sort by time
        [~, idx] = sort(raw(:,1)); 
        raw = raw(idx,:);
        
        % ================================================
        % Basic columns (1-12) - WITH UNIT CONVERSION
        % ================================================
        b.name = sheets{i}; 
        b.t = raw(:,1);              % Time (h)
        b.x = raw(:,2);              % Xv (10^6 cells/mL) - NO conversion
        b.via = raw(:,3)/100;        % Viability (convert % to fraction)
        
        % CONVERT g/L to mM
        b.gln = raw(:,4) / MW_GLN * 1000;  % GLN: g/L → mM
        b.glt = raw(:,5);                   % GLT (g/L) - keep as is
        b.glc = raw(:,6) / MW_GLC * 1000;  % GLC: g/L → mM
        b.lac = raw(:,7) / MW_LAC * 1000;  % LAC: g/L → mM
        b.nh4 = raw(:,8) / MW_NH4 * 1000;  % NH4: g/L → mM
        
        % Columns 9-11: ignored
        b.vol = raw(:,12) * 1000;           % Volume: L → mL (CONVERTED!)
        
        % ================================================
        % Feeding columns (13-18)
        % ================================================
        % Column 13: Pure glucose feed
        if size(raw,2) >= 13
            b.feed_glc_vol = raw(:,13) * 1000;  % L → mL
        else
            b.feed_glc_vol = zeros(size(b.t));
            warning('Column 13 (Glucose feed) not found. Using zeros.');
        end
        
        % Columns 14-18: Feed media 1-5 (variable number)
        feed_cols = 14:18;
        feed_data = [];
        n_feeds = 0;
        
        for col = feed_cols
            if size(raw,2) >= col
                feed_vol = raw(:,col) * 1000;  % L → mL
                % Check if this feed is used (has non-zero values)
                if any(feed_vol > 0, 'all') || any(~isnan(feed_vol), 'all')
                    n_feeds = n_feeds + 1;
                    feed_data(:, n_feeds) = feed_vol;
                end
            end
        end
        
        if n_feeds == 0
            feed_data = zeros(size(b.t), 1);
            n_feeds = 1;
            warning('No feed media columns found. Creating dummy feed.');
        end
        
        b.n_feeds = n_feeds;
        b.feed_vols = feed_data;  % [n_time x n_feeds]
        
        % Replace NaN with 0
        b.feed_glc_vol(isnan(b.feed_glc_vol)) = 0;
        b.feed_vols(isnan(b.feed_vols)) = 0;
        
        % ================================================
        % Feed compositions from Main parameter (ALREADY IN mM)
        % ================================================
        b.feed_glc = feed_glc_all(1:n_feeds);  % [n_feeds x 1] (mM)
        b.feed_gln = feed_gln_all(1:n_feeds);  % [n_feeds x 1] (mM)
        b.glc_stock = glc_stock;               % mM
        
        fprintf('   ✓ %d timepoints, %d feed media detected\n', length(b.t), n_feeds);
        fprintf('   ✓ Converted to mM: GLC, GLN, LAC, NH4\n');
        
        batches{i} = b;
    end
    
    fprintf('\n✅ Successfully loaded %d batches with unit conversion (g/L → mM)\n', nB);
end

function x = toNum(v)
    % Convert cell or string to number
    if isnumeric(v)
        x = v;
    else
        x = str2double(v); 
        if isnan(x)
            x = 0; 
        end
    end
end
