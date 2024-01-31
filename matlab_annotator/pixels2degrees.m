function [alphaH,alphaV] = pixels2degrees(pixels,d,pixelDimensions,stimulusDimensions)
% Converts pixels to degrees
% Input:  pixels (vector or scalar)
% Output: degrees of visual angle
% Observer that this conversion is not approariate close the the edges of
% a (wide) screen.
%--------------------------------------------------------------------------

meterPerPixelH = pixels.*(stimulusDimensions(1)/pixelDimensions(1));
meterPerPixelV = pixels.*(stimulusDimensions(2)/pixelDimensions(2));

alphaH = 180/pi*(2*atan(meterPerPixelH/d/2));
alphaV = 180/pi*(2*atan(meterPerPixelV/d/2));

% meterPerPixel = pixels*stimulusDimensions./pixelDimensions;
% alphaH = 180/pi*(2*atan(meterPerPixel(1)/d/2));
% alphaV = 180/pi*(2*atan(meterPerPixel(2)/d/2));

