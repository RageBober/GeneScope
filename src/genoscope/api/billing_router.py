"""
Billing API endpoints
Handles subscriptions, payments, and usage tracking
"""

from fastapi import APIRouter, HTTPException, Depends, Request, Header
from fastapi.responses import JSONResponse
from typing import Optional, Dict, Any
from datetime import datetime
import os
import logging

from ..billing.stripe_integration import (
    StripeManager,
    SubscriptionPlans,
    UsageTracker,
    BillingWebhooks
)

logger = logging.getLogger(__name__)

# Create router
billing_router = APIRouter(prefix="/api/billing", tags=["billing"])

# Initialize services
stripe_manager = StripeManager()
usage_tracker = UsageTracker()

# Webhook secret from Stripe dashboard
WEBHOOK_SECRET = os.getenv("STRIPE_WEBHOOK_SECRET", "whsec_test_...")

@billing_router.get("/plans")
async def get_subscription_plans():
    """Get all available subscription plans"""
    plans = SubscriptionPlans.get_all_plans()
    
    # Format for frontend
    formatted_plans = []
    for plan in plans:
        formatted_plans.append({
            "id": plan["id"],
            "name": plan["name"],
            "price": plan["price"],
            "price_display": f"${plan['price']}/month" if isinstance(plan["price"], int) else "Contact us",
            "limits": plan["limits"],
            "recommended": plan["id"] == "professional",
            "features": {
                "analyses_per_month": plan["limits"]["analyses_per_month"],
                "storage": f"{plan['limits']['storage_gb']} GB" if isinstance(plan["limits"]["storage_gb"], int) else "Unlimited",
                "api_calls": f"{plan['limits']['api_calls']:,}" if isinstance(plan["limits"]["api_calls"], int) else "Unlimited",
                "max_file_size": f"{plan['limits']['vcf_size_mb']} MB" if isinstance(plan["limits"]["vcf_size_mb"], int) else "Unlimited",
                "features_list": plan["limits"]["features"]
            }
        })
    
    return {"plans": formatted_plans}

@billing_router.post("/subscribe")
async def create_subscription(
    plan_id: str,
    user_email: str,
    user_name: str
):
    """Create a new subscription"""
    try:
        # Get plan details
        plan = SubscriptionPlans.get_plan(plan_id)
        if not plan:
            raise HTTPException(status_code=404, detail="Plan not found")
        
        if plan["id"] == "free":
            # Free plan doesn't need Stripe
            return {
                "status": "success",
                "message": "Free plan activated",
                "plan": plan["name"]
            }
        
        if plan["id"] == "enterprise":
            # Enterprise needs custom handling
            return {
                "status": "contact_sales",
                "message": "Please contact sales for enterprise pricing",
                "contact_email": "sales@genoscope.com"
            }
        
        # Create Stripe customer
        customer_id = stripe_manager.create_customer(
            email=user_email,
            name=user_name,
            metadata={"plan": plan_id}
        )
        
        # Create subscription
        subscription = stripe_manager.create_subscription(
            customer_id=customer_id,
            price_id=plan["stripe_price_id"],
            trial_days=14
        )
        
        return {
            "status": "success",
            "subscription": subscription,
            "message": f"Subscription to {plan['name']} plan created",
            "trial_ends": subscription.get("trial_end")
        }
        
    except Exception as e:
        logger.error(f"Subscription creation failed: {e}")
        raise HTTPException(status_code=400, detail=str(e))

@billing_router.post("/cancel-subscription")
async def cancel_subscription(
    subscription_id: str,
    immediate: bool = False
):
    """Cancel a subscription"""
    try:
        result = stripe_manager.cancel_subscription(subscription_id, immediate)
        
        return {
            "status": "success",
            "cancellation": result,
            "message": "Subscription cancelled successfully"
        }
        
    except Exception as e:
        logger.error(f"Subscription cancellation failed: {e}")
        raise HTTPException(status_code=400, detail=str(e))

@billing_router.post("/update-subscription")
async def update_subscription(
    subscription_id: str,
    new_plan_id: str
):
    """Update subscription to a different plan"""
    try:
        # Get new plan details
        new_plan = SubscriptionPlans.get_plan(new_plan_id)
        if not new_plan or not new_plan.get("stripe_price_id"):
            raise HTTPException(status_code=404, detail="Invalid plan")
        
        result = stripe_manager.update_subscription(
            subscription_id=subscription_id,
            new_price_id=new_plan["stripe_price_id"]
        )
        
        return {
            "status": "success",
            "update": result,
            "message": f"Subscription updated to {new_plan['name']} plan"
        }
        
    except Exception as e:
        logger.error(f"Subscription update failed: {e}")
        raise HTTPException(status_code=400, detail=str(e))

@billing_router.get("/usage")
async def get_usage_summary(user_id: str):
    """Get current month usage summary"""
    try:
        usage = usage_tracker.get_usage_summary(user_id)
        
        # Get user's current plan (mock for now)
        current_plan = SubscriptionPlans.PROFESSIONAL
        
        # Check limits
        limits_status = usage_tracker.check_limits(user_id, current_plan)
        
        return {
            "usage": usage,
            "limits": limits_status,
            "plan": current_plan["name"],
            "billing_period": usage["month"]
        }
        
    except Exception as e:
        logger.error(f"Failed to get usage summary: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@billing_router.post("/track-usage")
async def track_usage(
    user_id: str,
    analysis_type: str,
    size_mb: float
):
    """Track usage for billing purposes"""
    try:
        usage_update = usage_tracker.track_analysis(
            user_id=user_id,
            analysis_type=analysis_type,
            size_mb=size_mb
        )
        
        return {
            "status": "tracked",
            "current_usage": usage_update
        }
        
    except Exception as e:
        logger.error(f"Failed to track usage: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@billing_router.post("/create-payment")
async def create_one_time_payment(
    amount: int,
    description: str
):
    """Create a one-time payment for additional services"""
    try:
        payment = stripe_manager.create_payment_intent(
            amount=amount,
            currency="usd"
        )
        
        return {
            "status": "created",
            "payment": payment,
            "description": description
        }
        
    except Exception as e:
        logger.error(f"Failed to create payment: {e}")
        raise HTTPException(status_code=400, detail=str(e))

@billing_router.get("/portal-url")
async def get_customer_portal(
    customer_id: str,
    return_url: str = "http://localhost:3000/billing"
):
    """Get Stripe customer portal URL"""
    try:
        portal_url = stripe_manager.get_customer_portal_url(
            customer_id=customer_id,
            return_url=return_url
        )
        
        return {
            "portal_url": portal_url,
            "message": "Redirecting to billing portal..."
        }
        
    except Exception as e:
        logger.error(f"Failed to create portal session: {e}")
        raise HTTPException(status_code=400, detail=str(e))

@billing_router.post("/webhook")
async def stripe_webhook(
    request: Request,
    stripe_signature: str = Header(None)
):
    """Handle Stripe webhooks"""
    try:
        payload = await request.body()
        
        result = BillingWebhooks.handle_webhook(
            payload=payload,
            signature=stripe_signature,
            webhook_secret=WEBHOOK_SECRET
        )
        
        return JSONResponse(content=result, status_code=200)
        
    except ValueError as e:
        logger.error(f"Invalid webhook payload: {e}")
        raise HTTPException(status_code=400, detail="Invalid payload")
    except Exception as e:
        logger.error(f"Webhook processing failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@billing_router.get("/check-limits")
async def check_user_limits(
    user_id: str,
    action: str  # "analysis", "storage", "api_call"
):
    """Check if user can perform an action based on plan limits"""
    try:
        # Get user's plan (mock for now)
        user_plan = SubscriptionPlans.PROFESSIONAL
        
        limits = usage_tracker.check_limits(user_id, user_plan)
        
        can_proceed = {
            "analysis": limits["can_analyze"],
            "storage": limits["can_store"],
            "api_call": limits["can_call_api"]
        }.get(action, False)
        
        return {
            "allowed": can_proceed,
            "limits": limits["usage"],
            "plan": user_plan["name"],
            "upgrade_url": "/billing/upgrade" if not can_proceed else None
        }
        
    except Exception as e:
        logger.error(f"Failed to check limits: {e}")
        raise HTTPException(status_code=500, detail=str(e))
