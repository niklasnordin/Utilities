;; Make control-help copy a region on a Mac keyboard
;; This is to match control-instert under Linux
(global-set-key [(control help)] 'copy-region-as-kill)
;; Make shift-help paste the clipboard on a Mac keyboard
;; This is to match shift-instert under Linux
(global-set-key [(shift help)] 'yank-clipboard-selection)

(add-hook 'c-mode-common-hook 'openfoam-hgw-c-mode-hook)
 
;;(namespace-open  . 0)
;;(namespace-close . 0)
;;(innamespace     . 0)
