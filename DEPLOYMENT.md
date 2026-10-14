# Website deployment

Changes pushed to `main` deploy to the `dia-health-one` Cloudflare Pages project in the shared **Dia Health One** account. The Pages project's production branch must be configured as `main`. Other pushed branches create preview deployments. Pull requests validate website packaging without receiving deployment credentials.

The GitHub Actions workflow uploads `index.html`, optional `assets/` and `images/` directories, and optional `favicon.ico`, `robots.txt`, and `sitemap.xml`. Repository documentation and workflow files are excluded. Add other public files explicitly to the workflow when needed.

GitHub stores `CLOUDFLARE_API_TOKEN` as an encrypted Actions secret. The account-owned token grants Pages Write access only to account `009bb104fa671722398b10841f018688`. Either company Super Administrator can revoke or replace it through Cloudflare's Account API tokens page. Set its replacement with `gh secret set CLOUDFLARE_API_TOKEN --repo pallagip/dia-health-one`; enter the value through stdin, never commit it.

Deployments use a pinned Wrangler version and pinned GitHub actions. CI reads the unique deployment URL over HTTPS and compares its content with the committed page. The Pages production hostname is `dia-health-one.pages.dev`; custom-domain activation is a separate DNS/SSL check.

Domain registration and email remain at GoDaddy. Keep MX, SPF, DKIM, DMARC and autodiscover records intact and DNS-only. Do not enable Cloudflare Email Routing while GoDaddy provides the mailbox.

The original GitHub Pages workflow is retained during migration as a rollback option. To roll back website hosting after DNS migration, restore GitHub Pages' four A records, four AAAA records, and `www` CNAME. DNS backup and migration receipts are maintained outside this public website repository.
