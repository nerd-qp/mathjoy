;;(eval . (concat (file-name-as-directory default-directory) "build"))
;;((nil . ((helm-make-build-dir . "/home/king-kong/my-repo/hello-world/homework/偏微分方程数值解法/build"))))
;; ((nil . ((helm-make-build-dir . "build/")
;;          (helm-make-arguments . "-j8")
;;          (projectile-project-root . "/home/king-kong/my-repo/hello-world/homework/偏微分方程数值解法/"))))
;; (projectile-project-root . "/home/king-kong/my-repo/hello-world/homework/偏微分方程数值解法/code/")

((nil
  (helm-make-build-dir . "build")
  (helm-make-arguments . "-j8"))
 (c++-mode
  (flycheck-gcc-language-standard . "c++11")))

;;((nil . ((eval . (setq helm-make-build-dir (concat (file-name-as-directory default-directory) "build")))))
