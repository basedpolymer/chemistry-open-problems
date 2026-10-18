-- Cibles des renvois vers les figures, sans changer images, légendes ou placement.
-- Stocké avec le gabarit pour fonctionner aussi dans le dossier de partage.
function Pandoc(document)
  if not document.meta['cr-structured'] then return document end
  return document:walk({
    Figure = function(figure)
      if figure.identifier ~= '' then return figure end
      local source
      figure:walk({Image = function(image)
        source = source or image.src
      end})
      if source then
        local stem = source:match('([^/\\]+)%.[^%.]+$')
        if stem then figure.identifier = 'fig:' .. stem end
      end
      return figure
    end
  })
end
