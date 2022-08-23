# Administering the Lab Google Cloud Project

## Adding Users

Use the Google Cloud console's [IAM & Admin / IAM](https://console.cloud.google.com/iam-admin/iam?project=allen-discovery-center-mcovert) page to `+ADD` a new user with the below Roles (permissions).

**NOTE:** More narrow permissions would be more secure, but it'd take time and experiments to figure out what's really needed. E.g. `Logging / Logs Viewer` would probably be sufficient for debugging use, but then the person would need additional permissions to clean up old logs and maybe do some configuration.

User permissions/roles:
* Basic / Editor
* Cloud IAP / IAP-secured Tunnel User
* Logging / Logging Admin

Admin user additional permissions:
* Cloud IAP / IAP Policy Admin
* Cloud IAP / IAP Settings Admin
* Owner
