function gradientColors = createColorGradient(colorNames, numColors)
    % Get RGB values from the color library
    colors = cell(length(colorNames), 1);
    for i = 1:length(colorNames)
        colors{i} = DcolorLibrary(colorNames{i});
    end
    
    % Create a gradient between colors
    gradientColors = zeros(numColors, 3);
    for i = 1:3 % For R, G, B channels
        channel = linspace(0, 1, length(colors));
        channelGradient = zeros(numColors, 1);
        for j = 1:numColors
            pos = (j-1)/(numColors-1) * (length(colors) - 1) + 1;
            intPos = floor(pos);
            fracPos = pos - intPos;
            if intPos < length(colors)
                channelGradient(j) = (1 - fracPos) * colors{intPos}(i) + fracPos * colors{intPos + 1}(i);
            else
                channelGradient(j) = colors{end}(i);
            end
        end
        gradientColors(:, i) = channelGradient;
    end
end
