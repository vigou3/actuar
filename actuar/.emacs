;;;  -*- emacs-lisp -*-

;;;
;;; Fichier d'initialisation de GNU Emacs de Vincent Goulet
;;; (Version Windows simplifée)
;;;


;;; ==================================
;;;  Fichier de configuration externe
;;; ==================================
;(setq custom-file
;      (expand-file-name ".emacs-custom" "~"))
;(load-file custom-file)


;;; ============
;;;  LaTeX-mode
;;; ============
(setq-default TeX-auto-parse-length 200)
(setq-default TeX-master nil)
(setq LaTeX-default-options "12pt")

(add-hook 'LaTeX-mode-hook 'my-latex-mode)

(defun my-latex-mode ()
  (latex-math-mode)
  (turn-on-reftex)
  ;(setq ispell-extra-args '("-t"))
  (define-key LaTeX-math-keymap [?` (?_)] 'LaTeX-math-bar))


;;; =====
;;;  ESS
;;; =====
(setq-default inferior-R-program-name 
              "c:/documents and settings/lppou/r/r-2.3.0/bin/rterm.exe")

(require 'ess-site)

(defun Rnw-mode ()
  (require 'ess-noweb)
  (noweb-mode)
  (if (fboundp 'R-mode)
      (setq noweb-default-code-mode 'R-mode)))

(add-hook 'ess-mode-hook 'my-ess-options)
(add-hook 'inferior-ess-mode-hook 'my-iess-keybindings)

(defun my-ess-options ()
  (ess-set-style 'C++)
  (column-number-mode t)
  (define-key ess-mode-map [(meta backspace)] 'backward-kill-word)
  (define-key ess-mode-map [(control ?c) (?;)] 'comment-region)
  (define-key ess-mode-map [(control ?c) (?:)] 'uncomment-region))
  
(defun my-iess-keybindings ()
  (define-key inferior-ess-mode-map [(control ?a)] 'comint-bol)
  (define-key inferior-ess-mode-map [home] 'comint-bol))


;;; ==========
;;;  cc-mode
;;; ==========
(add-hook 'c-mode-common-hook 'my-cc-mode)

(defun my-cc-mode ()
  (font-lock-mode)
  (font-lock-fontify-buffer)
  (auto-fill-mode)
  (c-toggle-auto-state)
  (c-toggle-hungry-state)
  (c-set-style "bsd")
  (setq c-basic-offset 4))


;;; ===========
;;;  text-mode
;;; ===========
(add-hook 'text-mode-hook 'turn-on-auto-fill)


;;; ========
;;;  Ispell
;;; ========
;(require 'ispell)
;(setq ispell-dictionary "francais")


;;; ==========================================================
;;;  Nouvelles fonctions non rattachées à un mode particulier
;;; ==========================================================

;; =================================
;; Fonction pour dupliquer une ligne
;; =================================
(defun dup-line ()
  "Duplicates the line on which point lies."
  (interactive)
  (save-excursion
    (beginning-of-line)
    (let ((begin (point)))
      (forward-line)
      (copy-region-as-kill begin (point))
      (yank)
      (forward-line -1)
      (back-to-indentation))))

;;; ===============
;;;  Look and feel 
;;; ===============
(blink-cursor-mode nil)                  ; curseur ne clignote pas
(setq-default cursor-type 'bar)          ; curseur étroit
(set-face-background 'cursor "#CC0000")  ; curseur rouge foncé

(global-font-lock-mode t)                ; colorisation du texte
(transient-mark-mode t)                  ; mode de sélection "normal"
(delete-selection-mode t)                ; entrée efface texte sélectionné
(setq-default mouse-yank-at-point t)     ; coller avec la souris
(show-paren-mode t)                      ; coupler les parenthèses
(setq-default case-fold-search t)        ; recherche sans égard à la casse

(setq default-major-mode 'text-mode)     ; mode par défaut

(set-language-environment "Latin-1")     ; langage avec accents
(add-untranslated-filesystem "~/actuar")

;;; ==================================
;;;  Polices de caractère et couleurs
;;; ==================================
(set-face-font 'default "-outline-Courier New-normal-r-normal-normal-*-140-96-96-c-120-iso8859-1")
(set-face-font 'bold "-outline-Courier New-bold-r-normal-normal-*-140-96-96-c-120-iso8859-1")
(set-face-font 'italic "-outline-Courier New-normal-i-normal-normal-*-140-96-96-c-120-iso8859-1")
(set-face-font 'bold-italic "-outline-Courier New-bold-i-normal-normal-*-140-96-96-c-120-iso8859-1")
(set-face-background 'mode-line "#EFEEE2")
(set-face-background 'region "lightskyblue1")
(set-face-foreground 'font-lock-function-name-face "Blue3")
(set-face-foreground 'font-lock-keyword-face "Orange")


;;; ==================
;;;  Auto-compression 
;;; ==================
(auto-compression-mode nil)
(auto-compression-mode t)


;;; ================
;;;  Iswitch buffer
;;; ================
(iswitchb-mode t)


;;; ==============================
;;;  Nouveaux keybindings globaux
;;; ==============================
(global-set-key [(meta ?z)] 'dup-line)
(global-set-key [f5] 'goto-line)
(global-set-key [f8] 'auto-fill-mode)
(global-set-key [f12] 'other-window)
(global-set-key [mouse-4] 'scroll-down)
(global-set-key [mouse-5] 'scroll-up)
(global-set-key [(meta delete)] 'kill-word)
